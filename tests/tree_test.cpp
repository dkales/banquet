#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../banquet.h"
#include "../tree.h"

TEST_CASE("Tree is constructed", "[tree]") {
  seed_t seed = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  banquet_salt_t salt = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  SeedTree tree(seed, 64, salt, 0);
  for (size_t idx = 0; idx < 64; idx++) {
    REQUIRE(tree.get_leaf(idx));
  }
}

TEST_CASE("Reveallist is constructed", "[tree]") {
  seed_t seed = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  banquet_salt_t salt = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  SeedTree tree(seed, 64, salt, 0);
  reveal_list_t reveal_list = tree.reveal_all_but(0);
  REQUIRE(reveal_list.second == 0);
  REQUIRE(reveal_list.first.size() == 6);
}

TEST_CASE("Reveallist can reconstruct tree", "[tree]") {
  seed_t seed = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  banquet_salt_t salt = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  SeedTree tree(seed, 64, salt, 0);
  reveal_list_t reveal_list = tree.reveal_all_but(0);
  SeedTree tree2(reveal_list, 64, salt, 0);

  REQUIRE(tree2.get_leaf(0).has_value() == false);
  for (size_t idx = 1; idx < 64; idx++) {
    REQUIRE(tree.get_leaf(idx).value() == tree2.get_leaf(idx).value());
  }
}