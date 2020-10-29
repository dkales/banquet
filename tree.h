#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <vector>

#include "banquet.h"

constexpr size_t SEED_SIZE = 16;
constexpr size_t DIGEST_SIZE = 32;

typedef std::array<uint8_t, SEED_SIZE> seed_t;
typedef std::array<uint8_t, DIGEST_SIZE> digest_t;

class SeedTree {
public:
  typedef std::pair<std::vector<seed_t>, size_t> reveal_list_t;

private:
  // data, layed out continously in memory root at [0], its two children at [1],
  // [2], in general node at [n], children at [2*n + 1], [2*n + 2]
  std::vector<std::optional<seed_t>> _data;
  std::vector<bool> _node_exists;
  size_t _num_leaves;
  size_t _num_total_nodes;

  std::pair<seed_t, seed_t> expandSeed(const seed_t &seed,
                                       const banquet_salt_t &salt,
                                       const size_t rep_idx,
                                       const size_t node_idx);
  inline bool node_exists(size_t idx) {
    if (idx >= _num_total_nodes)
      return false;

    return _node_exists[idx];
  };

public:
  // construct from given seed, expand into num_leaves small seeds
  SeedTree(seed_t seed, const size_t num_leaves, const banquet_salt_t &salt,
           const size_t rep_idx);
  // re-construct from reveallist, expand all known values
  SeedTree(const reveal_list_t &reveallist, const size_t num_leaves,
           const banquet_salt_t &salt, const size_t rep_idx);
  ~SeedTree() = default;

  reveal_list_t reveal_all_but(size_t leaf_idx);
  std::optional<seed_t> get_leaf(size_t leaf_idx);
};

class MerkleTree {
private:
  /* data */
public:
  typedef std::array<uint8_t, DIGEST_SIZE> Digest;

  // build a merkle tree from a list of digests
  MerkleTree(const std::vector<MerkleTree::Digest> &digests);
  ~MerkleTree();

  // get root digest
  Digest get_root();
};