#include "tree.h"
#include "macros.h"

#include <cassert>
extern "C"
{
#include "kdf_shake.h"
}

static size_t get_parent(size_t node)
{
    assert(node != 0);
    return ((node + 1) >> 1) - 1;
}

SeedTree::SeedTree(SeedTree::seed_t seed, size_t num_leaves, const banquet_salt_t &salt, const size_t rep_idx) : _data(), _node_exists(), _num_leaves(num_leaves)
{
    size_t tree_depth = 1 + ceil_log2(num_leaves);

    size_t num_total_nodes =
        ((1 << (tree_depth)) - 1) -
        ((1 << (tree_depth - 1)) - num_leaves); /* Num nodes in complete - number of missing leaves */

    _node_exists.resize(num_total_nodes);
    for (size_t i = num_total_nodes - num_leaves; i < num_total_nodes; i++)
    {
        _node_exists[i] = true;
    }
    auto exists = [num_total_nodes, this](size_t idx) -> bool {
        if (idx >= num_total_nodes)
            return false;

        return _node_exists[idx];
    };

    for (size_t i = num_total_nodes - num_leaves; i > 0; i--)
    {
        if (exists(2 * i + 1) || exists(2 * i + 2))
        {
            _node_exists[i] = true;
        }
    }
    _node_exists[0] = true;
    // push back root seed
    _data.resize(num_total_nodes);
    _data[0] = std::optional(seed);

    size_t last_non_leaf = get_parent(num_total_nodes - 1);
    for (size_t i = 0; i < last_non_leaf; i++)
    {
        auto [left, right] = expandSeed(_data[0].value(), salt, rep_idx, i);
        if (exists(2 * i + 1))
        {
            _data[2 * i + 1] = std::optional(left);
        }
        if (exists(2 * i + 2))
        {
            _data[2 * i + 2] = std::optional(left);
        }
    }
}

std::pair<SeedTree::seed_t, SeedTree::seed_t> SeedTree::expandSeed(const SeedTree::seed_t &seed, const banquet_salt_t &salt, const size_t rep_idx, const size_t node_idx)
{
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

SeedTree::SeedTree(SeedTree::reveal_list_t reveallist)
{
}

SeedTree::reveal_list_t SeedTree::reveal_all_but(size_t leaf_idx)
{
    return SeedTree::reveal_list_t();
}