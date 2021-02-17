#include "tree.h"
#include "macros.h"

#include <cassert>
extern "C" {
#include "kdf_shake.h"
}

namespace {

void expand_seed(const gsl::span<uint8_t> &seed, const banquet_salt_t &salt,
                 const size_t rep_idx, const size_t node_idx,
                 gsl::span<uint8_t> &out) {
  hash_context ctx;

  hash_init_prefix(&ctx, seed.size() * 2, HASH_PREFIX_1);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, rep_idx);
  hash_update_uint16_le(&ctx, node_idx);
  hash_final(&ctx);
  hash_squeeze(&ctx, out.data(), seed.size());
  if (out.size() == 2 * seed.size())
    hash_squeeze(&ctx, out.data() + seed.size(), seed.size());
  hash_clear(&ctx);
}

void expand_seed_x4(const gsl::span<uint8_t> &seed0,
                    const gsl::span<uint8_t> &seed1,
                    const gsl::span<uint8_t> &seed2,
                    const gsl::span<uint8_t> &seed3, const banquet_salt_t &salt,
                    const size_t rep_idx, const size_t node_idx,
                    gsl::span<uint8_t> &out0, gsl::span<uint8_t> &out1,
                    gsl::span<uint8_t> &out2, gsl::span<uint8_t> &out3) {
  std::array<uint8_t, 32> dummy;
  hash_context_x4 ctx;
  const size_t seed_size = seed0.size();

  hash_init_prefix_x4(&ctx, seed_size * 2, HASH_PREFIX_1);
  hash_update_x4_4(&ctx, seed0.data(), seed1.data(), seed2.data(), seed3.data(),
                   seed_size);
  hash_update_x4_1(&ctx, salt.data(), salt.size());
  hash_update_x4_uint16_le(&ctx, rep_idx);
  const uint16_t node_ids[4] = {(uint16_t)(node_idx), (uint16_t)(node_idx + 1),
                                (uint16_t)(node_idx + 2),
                                (uint16_t)(node_idx + 3)};
  hash_update_x4_uint16s_le(&ctx, node_ids);
  hash_final_x4(&ctx);
  hash_squeeze_x4_4(&ctx, out0.data(), out1.data(), out2.data(), out3.data(),
                    seed_size);
  uint8_t *outptr[4] = {dummy.data(), dummy.data(), dummy.data(), dummy.data()};
  if (out0.size() == 2 * seed_size)
    outptr[0] = out0.data() + seed_size;
  if (out1.size() == 2 * seed_size)
    outptr[1] = out1.data() + seed_size;
  if (out2.size() == 2 * seed_size)
    outptr[2] = out2.data() + seed_size;
  if (out3.size() == 2 * seed_size)
    outptr[3] = out3.data() + seed_size;
  hash_squeeze_x4(&ctx, outptr, seed_size);
  hash_clear(&ctx);
}
size_t get_parent(size_t node) {
  assert(node != 0);
  return ((node + 1) >> 1) - 1;
}
} // namespace

SeedTree::SeedTree(const std::vector<uint8_t> &seed, const size_t num_leaves,
                   const banquet_salt_t &salt, const size_t rep_idx)
    : _data(), _node_exists(), _node_has_value(), _num_leaves(num_leaves) {
  size_t tree_depth = 1 + ceil_log2(num_leaves);

  _num_total_nodes =
      ((1 << (tree_depth)) - 1) -
      ((1 << (tree_depth - 1)) -
       num_leaves); /* Num nodes in complete - number of missing leaves */

  _node_exists.resize(_num_total_nodes);
  _node_has_value.resize(_num_total_nodes);
  for (size_t i = _num_total_nodes - num_leaves; i < _num_total_nodes; i++) {
    _node_exists[i] = true;
  }

  for (size_t i = _num_total_nodes - num_leaves; i > 0; i--) {
    if (node_exists(2 * i + 1) || node_exists(2 * i + 2)) {
      _node_exists[i] = true;
    }
  }
  _node_exists[0] = true;
  _node_has_value[0] = true;
  // push back root seed
  _seed_size = seed.size();
  _data.resize(_num_total_nodes * _seed_size);
  //_data[0] = seed;
  std::copy(std::begin(seed), std::end(seed), std::begin(_data));

  size_t last_non_leaf = get_parent(_num_total_nodes - 1);
  size_t i = 0;
  for (; i <= MIN(last_non_leaf, 2); i++) {
    if (!node_exists(i))
      continue;
    _node_has_value[2 * i + 1] = true;
    if (node_exists(2 * i + 2)) {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], 2 * _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
      _node_has_value[2 * i + 2] = true;
    } else {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
    }
  }
  // do some 4x iterations
  std::array<uint8_t, 64> dummy = {
      0,
  };
  for (; i < (last_non_leaf / 4) * 4; i += 4) {
    std::array<gsl::span<uint8_t>, 4> seeds = {
        gsl::span(dummy), gsl::span(dummy), gsl::span(dummy), gsl::span(dummy)};
    std::array<gsl::span<uint8_t>, 4> dsts = {
        gsl::span(dummy), gsl::span(dummy), gsl::span(dummy), gsl::span(dummy)};
    for (size_t j = 0; j < 4; j++) {
      if (node_exists(i + j)) {
        seeds[j] = gsl::span<uint8_t>(&_data[(i + j) * _seed_size], _seed_size);
        _node_has_value[2 * (i + j) + 1] = true;
        if (node_exists(2 * (i + j) + 2)) {
          dsts[j] =
              gsl::span(&_data[(2 * (i + j) + 1) * _seed_size], 2 * _seed_size);
          _node_has_value[2 * (i + j) + 2] = true;
        } else {
          dsts[j] =
              gsl::span(&_data[(2 * (i + j) + 1) * _seed_size], _seed_size);
        }
      }
    }
    expand_seed_x4(seeds[0], seeds[1], seeds[2], seeds[3], salt, rep_idx, i,
                   dsts[0], dsts[1], dsts[2], dsts[3]);
  }
  // do the remaining iterations
  for (; i <= last_non_leaf; i++) {
    if (!node_exists(i))
      continue;
    _node_has_value[2 * i + 1] = true;
    if (node_exists(2 * i + 2)) {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], 2 * _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
      _node_has_value[2 * i + 2] = true;
    } else {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
    }
  }
}

SeedTree::SeedTree(const reveal_list_t &reveallist, const size_t num_leaves,
                   const banquet_salt_t &salt, const size_t rep_idx)
    : _data(), _node_exists(), _node_has_value(), _num_leaves(num_leaves) {
  size_t tree_depth = 1 + ceil_log2(num_leaves);

  _num_total_nodes =
      ((1 << (tree_depth)) - 1) -
      ((1 << (tree_depth - 1)) -
       num_leaves); /* Num nodes in complete - number of missing leaves */

  _node_exists.resize(_num_total_nodes);
  _node_has_value.resize(_num_total_nodes);
  for (size_t i = _num_total_nodes - num_leaves; i < _num_total_nodes; i++) {
    _node_exists[i] = true;
  }

  for (size_t i = _num_total_nodes - num_leaves; i > 0; i--) {
    if (node_exists(2 * i + 1) || node_exists(2 * i + 2)) {
      _node_exists[i] = true;
    }
  }
  _node_exists[0] = true;
  _seed_size = reveallist.first.front().size();
  _data.resize(_num_total_nodes * _seed_size);

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
      path_idx++;
      continue;
    }
    size_t sibling = get_sibling(node);
    std::copy(std::begin(reveallist.first[path_idx]),
              std::end(reveallist.first[path_idx]),
              &_data[sibling * _seed_size]);
    _node_has_value[sibling] = true;
    path_idx++;
  }

  size_t last_non_leaf = get_parent(_num_total_nodes - 1);
  size_t i = 0;
  for (; i <= MIN(last_non_leaf, 2); i++) {
    if (!node_exists(i) || !node_has_value(i))
      continue;
    _node_has_value[2 * i + 1] = true;
    if (node_exists(2 * i + 2)) {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], 2 * _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
      _node_has_value[2 * i + 2] = true;
    } else {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
    }
  }
  // do some 4x iterations
  std::array<uint8_t, 64> dummy;
  std::fill(std::begin(dummy), std::end(dummy), 0);
  for (; i < (last_non_leaf / 4) * 4; i += 4) {
    std::array<gsl::span<uint8_t>, 4> seeds = {
        gsl::span(dummy.data(), _seed_size),
        gsl::span(dummy.data(), _seed_size),
        gsl::span(dummy.data(), _seed_size),
        gsl::span(dummy.data(), _seed_size)};
    std::array<gsl::span<uint8_t>, 4> dsts = {
        gsl::span(dummy.data(), 2 * _seed_size),
        gsl::span(dummy.data(), 2 * _seed_size),
        gsl::span(dummy.data(), 2 * _seed_size),
        gsl::span(dummy.data(), 2 * _seed_size)};
    for (size_t j = 0; j < 4; j++) {
      if (node_exists(i + j) && node_has_value(i + j)) {
        seeds[j] = gsl::span<uint8_t>(&_data[(i + j) * _seed_size], _seed_size);
        _node_has_value[2 * (i + j) + 1] = true;
        if (node_exists(2 * (i + j) + 2)) {
          dsts[j] =
              gsl::span(&_data[(2 * (i + j) + 1) * _seed_size], 2 * _seed_size);
          _node_has_value[2 * (i + j) + 2] = true;
        } else {
          dsts[j] =
              gsl::span(&_data[(2 * (i + j) + 1) * _seed_size], _seed_size);
        }
      }
    }
    expand_seed_x4(seeds[0], seeds[1], seeds[2], seeds[3], salt, rep_idx, i,
                   dsts[0], dsts[1], dsts[2], dsts[3]);
  }
  // do the remaining iterations
  for (; i <= last_non_leaf; i++) {
    if (!node_exists(i) || !node_has_value(i))
      continue;
    _node_has_value[2 * i + 1] = true;
    if (node_exists(2 * i + 2)) {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], 2 * _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
      _node_has_value[2 * i + 2] = true;
    } else {
      gsl::span dst(&_data[(2 * i + 1) * _seed_size], _seed_size);
      expand_seed(gsl::span<uint8_t>(&_data[i * _seed_size], _seed_size), salt,
                  rep_idx, i, dst);
    }
  }
}

reveal_list_t SeedTree::reveal_all_but(size_t leaf_idx) {
  // calculate path up to root for missing leaf
  std::vector<std::vector<uint8_t>> path;

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
      path.push_back(std::vector<uint8_t>(_seed_size));
      continue;
    }
    size_t sibling = get_sibling(node);
    std::vector<uint8_t> value;
    value.reserve(_seed_size);
    std::copy(&_data[sibling * _seed_size],
              &_data[sibling * _seed_size + _seed_size],
              std::back_inserter(value));
    path.push_back(value);
  }

  return std::make_pair(path, leaf_idx);
}

std::optional<gsl::span<uint8_t>> SeedTree::get_leaf(size_t leaf_idx) {
  assert(leaf_idx < _num_leaves);
  size_t real_leaf_idx = _num_total_nodes - _num_leaves + leaf_idx;
  if (!node_has_value(real_leaf_idx))
    return std::nullopt;
  return gsl::span<uint8_t>(&_data[real_leaf_idx * _seed_size], _seed_size);
}