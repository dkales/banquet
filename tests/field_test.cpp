#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../field.h"

TEST_CASE("Basic Arithmetic GF(2^32)", "[field]") {
  GF2_32 zero;
  GF2_32 one(1);
  GF2_32 x(2);
  GF2_32 x_2(4);

  REQUIRE(one + zero == one);
  REQUIRE(zero + one == one);
  REQUIRE(zero + zero == zero);
  REQUIRE(one + one == zero);

  REQUIRE(one * one == one);
  REQUIRE(zero * one == zero);
  REQUIRE(one * zero == zero);
  REQUIRE(zero * zero == zero);

  REQUIRE(x * one == x);
  REQUIRE(x * x == x_2);
}
TEST_CASE("Modular Arithmetic GF(2^32)", "[field]") {
  GF2_32 a, b;
  a.set_coeff(31);
  a.set_coeff(29);
  a.set_coeff(28);
  a.set_coeff(24);
  a.set_coeff(23);
  a.set_coeff(21);
  a.set_coeff(19);
  a.set_coeff(15);
  a.set_coeff(14);
  a.set_coeff(9);
  a.set_coeff(8);
  a.set_coeff(0);

  b.set_coeff(29);
  b.set_coeff(27);
  b.set_coeff(26);
  b.set_coeff(25);
  b.set_coeff(20);
  b.set_coeff(17);
  b.set_coeff(14);
  b.set_coeff(11);
  b.set_coeff(10);
  b.set_coeff(5);
  b.set_coeff(3);
  b.set_coeff(2);
  GF2_32 a_int(2980627201), b_int(772951084);
  REQUIRE(a == a_int);
  REQUIRE(b == b_int);

  GF2_32 ab(3895706975);
  GF2_32 a_2(2846443116);
  GF2_32 b_2(234130046);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(b * b == b_2);
}