#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../field.h"
using namespace field;

TEST_CASE("Basic Arithmetic in all fields", "[field]") {
  banquet_params_t params[] = {Banquet_L1_Param1, Banquet_L1_Param3,
                               Banquet_L1_Param4};
  for (auto param : params) {
    GF2E::init_extension_field(banquet_instance_get(param));
    GF2E zero;
    GF2E one(1);
    GF2E x(2);
    GF2E x_2(4);

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
}
TEST_CASE("Modular Arithmetic GF(2^32)", "[field]") {
  GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  GF2E a, b;
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
  GF2E a_int(2980627201), b_int(772951084);
  REQUIRE(a == a_int);
  REQUIRE(b == b_int);

  GF2E ab(3895706975);
  GF2E a_2(2846443116);
  GF2E b_2(234130046);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(b * b == b_2);
}

TEST_CASE("Modular Arithmetic GF(2^40)", "[field]") {
  GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param3));
  GF2E a(683080732428ULL), b(437065243549ULL);

  GF2E ab(916833885315ULL);
  GF2E a_2(224904587486ULL);
  GF2E b_2(153370336291ULL);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(b * b == b_2);
}

TEST_CASE("Modular Arithmetic GF(2^48)", "[field]") {
  GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param4));
  GF2E a(94834555057164ULL), b(119504161027501ULL);

  GF2E ab(130190305526807ULL);
  GF2E a_2(187414204277907ULL);
  GF2E b_2(119412372920018ULL);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(b * b == b_2);
}