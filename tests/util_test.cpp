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
