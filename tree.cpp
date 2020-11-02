#include "tree.h"
#include "macros.h"

#include <cassert>
extern "C" {
#include "kdf_shake.h"
}

static size_t get_parent(size_t node) {
  assert(node != 0);
  return ((node + 1) >> 1) - 1;
}

SeedTree::SeedTree(seed_t seed, const size_t num_leaves,
                   const banquet_salt_t &salt, const size_t rep_idx)
    : _data(), _node_exists(), _num_leaves(num_leaves) {
  size_t tree_depth = 1 + ceil_log2(num_leaves);

  _num_total_nodes =
      ((1 << (tree_depth)) - 1) -
      ((1 << (tree_depth - 1)) -
       num_leaves); /* Num nodes in complete - number of missing leaves */

  _node_exists.resize(_num_total_nodes);
  for (size_t i = _num_total_nodes - num_leaves; i < _num_total_nodes; i++) {
    _node_exists[i] = true;
  }

  for (size_t i = _num_total_nodes - num_leaves; i > 0; i--) {
    if (node_exists(2 * i + 1) || node_exists(2 * i + 2)) {
      _node_exists[i] = true;
    }
  }
  _node_exists[0] = true;
  // push back root seed
  _data.resize(_num_total_nodes);
  _data[0] = std::optional(seed);

  size_t last_non_leaf = get_parent(_num_total_nodes - 1);
  for (size_t i = 0; i <= last_non_leaf; i++) {
    auto [left, right] = expandSeed(_data[i].value(), salt, rep_idx, i);
    if (node_exists(2 * i + 1)) {
      _data[2 * i + 1] = std::optional(left);
    }
    if (node_exists(2 * i + 2)) {
      _data[2 * i + 2] = std::optional(right);
    }
  }
}

std::pair<seed_t, seed_t> SeedTree::expandSeed(const seed_t &seed,
                                               const banquet_salt_t &salt,
                                               const size_t rep_idx,
                                               const size_t node_idx) {
  hash_context ctx;
  std::array<uint8_t, SEED_SIZE> ret1;
  std::array<uint8_t, SEED_SIZE> ret2;

  hash_init_prefix(&ctx, DIGEST_SIZE, HASH_PREFIX_1);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, rep_idx);
  hash_update_uint16_le(&ctx, node_idx);
  hash_final(&ctx);
  hash_squeeze(&ctx, ret1.data(), SEED_SIZE);
  hash_squeeze(&ctx, ret2.data(), SEED_SIZE);
  hash_clear(&ctx);
  return std::make_pair(ret1, ret2);
}

SeedTree::SeedTree(const reveal_list_t &reveallist, const size_t num_leaves,
                   const banquet_salt_t &salt, const size_t rep_idx)
    : _data(), _node_exists(), _num_leaves(num_leaves) {
  size_t tree_depth = 1 + ceil_log2(num_leaves);

  _num_total_nodes =
      ((1 << (tree_depth)) - 1) -
      ((1 << (tree_depth - 1)) -
       num_leaves); /* Num nodes in complete - number of missing leaves */

  _node_exists.resize(_num_total_nodes);
  for (size_t i = _num_total_nodes - num_leaves; i < _num_total_nodes; i++) {
    _node_exists[i] = true;
  }

  for (size_t i = _num_total_nodes - num_leaves; i > 0; i--) {
    if (node_exists(2 * i + 1) || node_exists(2 * i + 2)) {
      _node_exists[i] = true;
    }
  }
  _node_exists[0] = true;
  _data.resize(_num_total_nodes);

  auto has_sibling = [this](size_t node) -> bool {
    if (!node_exists(node)) {
      return 0;
    }

    if ((node % 2 == 1) && !node_exists(node + 1)) {
      return 0;
    }

    return 1;
  };

  auto get_sibling = [this, has_sibling](size_t node) -> size_t {
    assert(node < _num_total_nodes);
    assert(node != 0);
    assert(has_sibling(node));
    if ((node % 2 == 1)) {
      if (node + 1 < _num_total_nodes) {
        return node + 1;
      } else {
        assert(!"getSibling: request for node with no sibling");
        return 0;
      }
    } else {
      return node - 1;
    }
  };

  size_t missing_leaf = reveallist.second;
  size_t first_leaf_idx = _num_total_nodes - _num_leaves;
  size_t path_idx = 0;
  for (size_t node = first_leaf_idx + missing_leaf; node != 0;
       node = get_parent(node)) {
    if (!has_sibling(node)) {
      continue;
    }
    size_t sibling = get_sibling(node);
    _data[sibling] = std::optional(reveallist.first[path_idx]);
    path_idx++;
  }

  size_t last_non_leaf = get_parent(_num_total_nodes - 1);
  for (size_t i = 0; i <= last_non_leaf; i++) {
    if (!_data[i].has_value())
      continue;
    auto [left, right] = expandSeed(_data[i].value(), salt, rep_idx, i);
    if (node_exists(2 * i + 1)) {
      _data[2 * i + 1] = std::optional(left);
    }
    if (node_exists(2 * i + 2)) {
      _data[2 * i + 2] = std::optional(right);
    }
  }
}

reveal_list_t SeedTree::reveal_all_but(size_t leaf_idx) {
  // calculate path up to root for missing leaf
  std::vector<seed_t> path;

  auto has_sibling = [this](size_t node) -> bool {
    if (!node_exists(node)) {
      return 0;
    }

    if ((node % 2 == 1) && !node_exists(node + 1)) {
      return 0;
    }

    return 1;
  };

  auto get_sibling = [this, has_sibling](size_t node) -> size_t {
    assert(node < _num_total_nodes);
    assert(node != 0);
    assert(has_sibling(node));
    if ((node % 2 == 1)) {
      if (node + 1 < _num_total_nodes) {
        return node + 1;
      } else {
        assert(!"getSibling: request for node with no sibling");
        return 0;
      }
    } else {
      return node - 1;
    }
  };
  size_t first_leaf_idx = _num_total_nodes - _num_leaves;
  for (size_t node = first_leaf_idx + leaf_idx; node != 0;
       node = get_parent(node)) {
    if (!has_sibling(node)) {
      continue;
    }
    size_t sibling = get_sibling(node);
    path.push_back(_data[sibling].value());
  }

  return std::make_pair(path, leaf_idx);
}

std::optional<seed_t> SeedTree::get_leaf(size_t leaf_idx) {
  assert(leaf_idx < _num_leaves);
  size_t real_leaf_idx = _num_total_nodes - _num_leaves + leaf_idx;
  return _data[real_leaf_idx];
}