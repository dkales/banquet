#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../banquet.h"
#include "../utils.h"
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>

TEST_CASE("Precomputed Lagrange Interpolation", "[util]") {
  const banquet_instance_t &instance = banquet_instance_get(Banquet_L1_Param1);
  utils::init_extension_field(instance);
  size_t dimension = 20;
  vec_GF2E x_values;
  vec_GF2E y_values1;
  vec_GF2E y_values2;
  for (size_t i = 0; i < dimension; i++) {
    x_values.append(random_GF2E());
    y_values1.append(random_GF2E());
    y_values2.append(random_GF2E());
  }

  // builtin interpolate
  GF2EX poly1 = interpolate(x_values, y_values1);
  GF2EX poly2 = interpolate(x_values, y_values2);
  // precomputed interpolate
  auto precomputation = utils::precompute_lagrange_polynomials(x_values);
  GF2EX poly1_with_precom =
      utils::interpolate_with_precomputation(precomputation, y_values1);
  GF2EX poly2_with_precom =
      utils::interpolate_with_precomputation(precomputation, y_values2);

  REQUIRE(poly1 == poly1_with_precom);
  REQUIRE(poly2 == poly2_with_precom);
}

TEST_CASE("Basic Lifting tests", "[util]") {
  const banquet_instance_t &instance = banquet_instance_get(Banquet_L1_Param1);
  utils::init_extension_field(instance);
  uint8_t a = 3;
  uint8_t b = 246;

  GF2E a_lifted, b_lifted;
  GF2X tmp;
  // a_lifted should be y^30 + y^23 + y^21 + y^18 + y^14 + y^13 + y^11 + y^9 +
  // y^7 + y^6 + y^5 + y^4 + y^3 + y + 1
  SetCoeff(tmp, 30);
  SetCoeff(tmp, 23);
  SetCoeff(tmp, 21);
  SetCoeff(tmp, 18);
  SetCoeff(tmp, 14);
  SetCoeff(tmp, 13);
  SetCoeff(tmp, 11);
  SetCoeff(tmp, 9);
  SetCoeff(tmp, 7);
  SetCoeff(tmp, 6);
  SetCoeff(tmp, 5);
  SetCoeff(tmp, 4);
  SetCoeff(tmp, 3);
  SetCoeff(tmp, 1);
  SetCoeff(tmp, 0);
  a_lifted = conv<GF2E>(tmp);
  clear(tmp);
  // b_lifted should be y^30 + y^29 + y^27 + y^26 + y^25 + y^19 + y^18 + y^17 +
  // y^14 + y^13 + y^12 + y^10 + y^8 + y^7 + y^4 + y^2 + y
  SetCoeff(tmp, 30);
  SetCoeff(tmp, 29);
  SetCoeff(tmp, 27);
  SetCoeff(tmp, 26);
  SetCoeff(tmp, 25);
  SetCoeff(tmp, 19);
  SetCoeff(tmp, 18);
  SetCoeff(tmp, 17);
  SetCoeff(tmp, 14);
  SetCoeff(tmp, 13);
  SetCoeff(tmp, 12);
  SetCoeff(tmp, 10);
  SetCoeff(tmp, 8);
  SetCoeff(tmp, 7);
  SetCoeff(tmp, 4);
  SetCoeff(tmp, 2);
  SetCoeff(tmp, 1);
  b_lifted = conv<GF2E>(tmp);
  GF2E one;
  set(one);

  REQUIRE(utils::lift_uint8_t(a) == a_lifted);
  REQUIRE(utils::lift_uint8_t(b) == b_lifted);
  REQUIRE(a_lifted * b_lifted == one);
}
