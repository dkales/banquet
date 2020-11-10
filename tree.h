#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <vector>

#include "types.h"

class SeedTree {
public:
private:
  // data, layed out continously in memory root at [0], its two children at [1],
  // [2], in general node at [n], children at [2*n + 1], [2*n + 2]
  std::vector<std::optional<std::vector<uint8_t>>> _data;
  std::vector<bool> _node_exists;
  size_t _num_leaves;
  size_t _num_total_nodes;

  std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
  expandSeed(const std::vector<uint8_t> &seed, const banquet_salt_t &salt,
             const size_t rep_idx, const size_t node_idx);
  inline bool node_exists(size_t idx) {
    if (idx >= _num_total_nodes)
      return false;

    return _node_exists[idx];
  };

public:
  // construct from given seed, expand into num_leaves small seeds
  SeedTree(std::vector<uint8_t> seed, const size_t num_leaves,
           const banquet_salt_t &salt, const size_t rep_idx);
  // re-construct from reveallist, expand all known values
  SeedTree(const reveal_list_t &reveallist, const size_t num_leaves,
           const banquet_salt_t &salt, const size_t rep_idx);
  ~SeedTree() = default;

  reveal_list_t reveal_all_but(size_t leaf_idx);
  std::optional<std::vector<uint8_t>> get_leaf(size_t leaf_idx);
};
