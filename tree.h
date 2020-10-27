#pragma once

#include <array>
#include <optional>
#include <vector>
#include <cstdint>
#include <cstdlib>

#include "banquet.h"

constexpr size_t SEED_SIZE = 16;
constexpr size_t DIGEST_SIZE = 32;

class SeedTree
{
    typedef std::array<uint8_t, SEED_SIZE> seed_t;
    typedef std::pair<std::vector<seed_t>, size_t> reveal_list_t;

private:
    // data, layed out continously in memory root at [0], its two children at [1], [2], in general node at [n], children at [2*n + 1], [2*n + 2]
    std::vector<std::optional<seed_t>> _data;
    std::vector<bool> _node_exists;
    size_t _num_leaves;

    std::pair<seed_t, seed_t> expandSeed(const seed_t &seed, const banquet_salt_t &salt, const size_t rep_idx, const size_t node_idx);

public:
    // construct from given seed, expand into num_leaves small seeds
    SeedTree(seed_t seed, size_t num_leaves, const banquet_salt_t &salt, const size_t rep_idx);
    // re-construct from reveallist, expand all known values
    SeedTree(reveal_list_t reveallist);
    ~SeedTree();

    reveal_list_t reveal_all_but(size_t leaf_idx);
    std::optional<seed_t> get_leaf(size_t leaf_idx);
};

class MerkleTree
{
private:
    /* data */
public:
    typedef std::array<uint8_t, DIGEST_SIZE> Digest;

    // build a merkle tree from a list of digests
    MerkleTree(const std::vector<MerkleTree::Digest> &digests);
    ~MerkleTree();

    //get root digest
    Digest get_root();
};