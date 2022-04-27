#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../field.h"
#include "utils.h"

#include <NTL/GF2EX.h>

TEST_CASE("Basic Arithmetic in all fields", "[field]") {
  banquet_params_t params[] = {Banquet_L1_Param1, Banquet_L1_Param3,
                               Banquet_L1_Param4};
  for (auto param : params) {
    field::GF2E::init_extension_field(banquet_instance_get(param));
    field::GF2E zero;
    field::GF2E one(1);
    field::GF2E x(2);
    field::GF2E x_2(4);

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
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  field::GF2E a, b;
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
  field::GF2E a_int(2980627201), b_int(772951084);
  REQUIRE(a == a_int);
  REQUIRE(b == b_int);

  field::GF2E ab(3895706975);
  field::GF2E a_2(2846443116);
  field::GF2E b_2(234130046);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(a.sqr() == a_2);
  REQUIRE(b * b == b_2);
  REQUIRE(b.sqr() == b_2);
}

TEST_CASE("Modular Arithmetic GF(2^40)", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param3));
  field::GF2E a(683080732428ULL), b(437065243549ULL);

  field::GF2E ab(916833885315ULL);
  field::GF2E a_2(224904587486ULL);
  field::GF2E b_2(153370336291ULL);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(a.sqr() == a_2);
  REQUIRE(b * b == b_2);
}

TEST_CASE("Modular Arithmetic GF(2^48)", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param4));
  field::GF2E a(94834555057164ULL), b(119504161027501ULL);

  field::GF2E ab(130190305526807ULL);
  field::GF2E a_2(187414204277907ULL);
  field::GF2E b_2(119412372920018ULL);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
  REQUIRE(a.sqr() == a_2);
  REQUIRE(b * b == b_2);
}

TEST_CASE("NTL to_bytes = custom to_bytes", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  utils::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  field::GF2E a, b;
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
  GF2X c, d;
  SetCoeff(c, 31);
  SetCoeff(c, 29);
  SetCoeff(c, 28);
  SetCoeff(c, 24);
  SetCoeff(c, 23);
  SetCoeff(c, 21);
  SetCoeff(c, 19);
  SetCoeff(c, 15);
  SetCoeff(c, 14);
  SetCoeff(c, 9);
  SetCoeff(c, 8);
  SetCoeff(c, 0);
  GF2E c_e = conv<GF2E>(c);

  SetCoeff(d, 29);
  SetCoeff(d, 27);
  SetCoeff(d, 26);
  SetCoeff(d, 25);
  SetCoeff(d, 20);
  SetCoeff(d, 17);
  SetCoeff(d, 14);
  SetCoeff(d, 11);
  SetCoeff(d, 10);
  SetCoeff(d, 5);
  SetCoeff(d, 3);
  SetCoeff(d, 2);
  GF2E d_e = conv<GF2E>(d);

  const GF2X &poly_rep_c = rep(c_e);
  std::vector<uint8_t> buffer_c(4);
  BytesFromGF2X(buffer_c.data(), poly_rep_c, buffer_c.size());
  const GF2X &poly_rep_d = rep(d_e);
  std::vector<uint8_t> buffer_d(4);
  BytesFromGF2X(buffer_d.data(), poly_rep_d, buffer_d.size());

  std::vector<uint8_t> buffer_a(4);
  a.to_bytes(buffer_a.data());
  std::vector<uint8_t> buffer_b(4);
  b.to_bytes(buffer_b.data());
  REQUIRE(buffer_a == buffer_c);
  REQUIRE(buffer_b == buffer_d);
}

TEST_CASE("NTL to custom conversion", "[field]") {
  banquet_params_t params[] = {Banquet_L1_Param1, Banquet_L1_Param3,
                               Banquet_L1_Param4};
  for (auto param : params) {
    utils::init_extension_field(banquet_instance_get(param));
    field::GF2E::init_extension_field(banquet_instance_get(param));
    field::GF2E a, b;
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

    field::GF2E ab = a * b;
    GF2E a_ntl = utils::custom_to_ntl(a);
    GF2E ab_ntl = utils::custom_to_ntl(ab);
    GF2E b_ntl = ab_ntl / a_ntl;
    field::GF2E b2 = utils::ntl_to_custom(b_ntl);
    REQUIRE(b == b2);
  }
}

TEST_CASE("NTL inverse == custom", "[field]") {
  utils::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  field::GF2E a;
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

  field::GF2E b = a.inverse();
  field::GF2E c = utils::ntl_to_custom(inv(utils::custom_to_ntl(a)));
  REQUIRE(b == c);
  REQUIRE(a * b == field::GF2E(1));
}

TEST_CASE("Constant time inverse == custom", "[field]") {

  // Checking for GF2_32
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  field::GF2E b;
  b.set_coeff(31);
  REQUIRE(b.inverse_const_time() == b.inverse());

  // Checking for GF2_40
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param3));
  field::GF2E c;
  c.set_coeff(39);
  REQUIRE(c.inverse_const_time() == c.inverse());

  // Checking for GF2_48
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param4));
  field::GF2E d;
  d.set_coeff(47);
  REQUIRE(d.inverse_const_time() == d.inverse());
}

TEST_CASE("RANDOM TESTS") {

  // field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  // field::GF2E a;
  // a.set_coeff(38);
  // a.set_coeff(31);

  // field::GF2E a1 = a.sqr();
  // field::GF2E a2 = a * a;

  // uint64_t round = 10000000;

  // __m128i a = _mm_set_epi64x(0x00, 0xbbbbbbbbbbbbbbbb);

  // REQUIRE(field::reduce_clmul_GF2(a) == field::reduce_GF2_barret(a));
  // REQUIRE(field::reduce_clmul_GF2(a) == field::reduce_GF2(a));

  // for (uint64_t i = 0; i < 1000; ++i) {
  //   field::reduce_clmul_GF2(a);
  // }
  // auto start = std::chrono::system_clock::now();
  // for (uint64_t i = 0; i < round; ++i) {
  //   field::reduce_clmul_GF2(a);
  // }
  // auto end = std::chrono::system_clock::now() - start;
  // std::cout << "clmul - " << field::reduce_clmul_GF2(a) << " - ";
  // std::cout << std::dec << end / std::chrono::milliseconds(1);
  // std::cout << "ms" << std::endl;

  // for (uint64_t i = 0; i < 1000; ++i) {
  //   field::reduce_GF2_barret(a);
  // }
  // start = std::chrono::system_clock::now();
  // for (uint64_t i = 0; i < round; ++i) {
  //   field::reduce_GF2_barret(a);
  // }
  // end = std::chrono::system_clock::now() - start;
  // std::cout << "barret - " << field::reduce_GF2_barret(a) << " - ";
  // std::cout << std::dec << end / std::chrono::milliseconds(1);
  // std::cout << "mms" << std::endl;

  // for (uint64_t i = 0; i < 1000; ++i) {
  //   field::reduce_GF2(a);
  // }
  // start = std::chrono::system_clock::now();
  // for (uint64_t i = 0; i < round; ++i) {
  //   field::reduce_GF2(a);
  // }
  // end = std::chrono::system_clock::now() - start;
  // std::cout << "normal - " << field::reduce_GF2(a) << " - ";
  // std::cout << std::dec << end / std::chrono::milliseconds(1);
  // std::cout << "ms" << std::endl;
}

TEST_CASE("NTL interpolation == custom", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  utils::init_extension_field(banquet_instance_get(Banquet_L1_Param1));

  std::vector<field::GF2E> a = field::get_first_n_field_elements(100);
  vec_GF2E b = utils::get_first_n_field_elements(100);
  for (size_t i = 0; i < 100; i++) {
    REQUIRE(a[i] == utils::ntl_to_custom(b[i]));
  }
  std::vector<field::GF2E> a_from_roots = field::build_from_roots(a);
  GF2EX b_from_roots = BuildFromRoots(b);
  REQUIRE(a_from_roots.size() == (size_t)b_from_roots.rep.length());
  for (size_t j = 0; j < a_from_roots.size(); j++) {
    REQUIRE(a_from_roots[j] == utils::ntl_to_custom(b_from_roots[j]));
  }

  std::vector<std::vector<field::GF2E>> a_lag =
      field::precompute_lagrange_polynomials_slow(a);
  std::vector<GF2EX> b_lag = utils::precompute_lagrange_polynomials(b);

  REQUIRE(a_lag.size() == b_lag.size());
  for (size_t i = 0; i < a_lag.size(); i++) {
    REQUIRE(a_lag[i].size() == (size_t)b_lag[i].rep.length());
    for (size_t j = 0; j < a_lag[i].size(); j++) {
      REQUIRE(a_lag[i][j] == utils::ntl_to_custom(b_lag[i][j]));
    }
  }
}

TEST_CASE("optmized custom == custom interpolation", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  const size_t ROOT_SIZE = 128;

  std::vector<field::GF2E> x = field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<field::GF2E> y = field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<std::vector<field::GF2E>> a_lag =
      field::precompute_lagrange_polynomials_slow(x);
  std::vector<field::GF2E> result =
      field::interpolate_with_precomputation(a_lag, y);

  std::vector<field::GF2E> x_opti =
      field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<field::GF2E> y_opti =
      field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<std::vector<field::GF2E>> x_lag =
      field::precompute_lagrange_polynomials(x_opti);
  std::vector<field::GF2E> result_optim =
      field::interpolate_with_precomputation(x_lag, y_opti);

  REQUIRE(result == result_optim);
}

TEST_CASE("fast interpolation == optmized custom interpolation", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param1));
  const size_t ROOT_SIZE = 128;

  std::vector<field::GF2E> x_opti =
      field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<field::GF2E> y_opti =
      field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<std::vector<field::GF2E>> x_lag =
      field::precompute_lagrange_polynomials(x_opti);

  std::vector<field::GF2E> result_optim =
      field::interpolate_with_precomputation(x_lag, y_opti);

  std::vector<field::GF2E> x_fast =
      field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<field::GF2E> y_fast =
      field::get_first_n_field_elements(ROOT_SIZE);
  std::vector<field::GF2E> precomputed_denominator =
      field::precompute_denominator(x_fast);
  std::vector<std::vector<field::GF2E>> precomputed_x_minus_xi;
  field::set_x_minus_xi_poly_size(precomputed_x_minus_xi, x_fast.size());
  field::precompute_x_minus_xi_poly_splits(x_fast, precomputed_x_minus_xi);

  std::vector<field::GF2E> result_fast = field::interpolate_with_recurrsion(
      y_fast, precomputed_denominator, precomputed_x_minus_xi, 0, x_fast.size(),
      0, precomputed_x_minus_xi.size());

  REQUIRE(result_fast == result_optim);
}
