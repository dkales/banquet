#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../field.h"
#include "../utils.h"

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
  REQUIRE(b * b == b_2);
}

TEST_CASE("Modular Arithmetic GF(2^40)", "[field]") {
  field::GF2E::init_extension_field(banquet_instance_get(Banquet_L1_Param3));
  field::GF2E a(683080732428ULL), b(437065243549ULL);

  field::GF2E ab(916833885315ULL);
  field::GF2E a_2(224904587486ULL);
  field::GF2E b_2(153370336291ULL);
  REQUIRE(a * b == ab);
  REQUIRE(a * a == a_2);
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
  std::vector<uint8_t> buffer_c(8);
  BytesFromGF2X(buffer_c.data(), poly_rep_c, buffer_c.size());
  const GF2X &poly_rep_d = rep(d_e);
  std::vector<uint8_t> buffer_d(8);
  BytesFromGF2X(buffer_d.data(), poly_rep_d, buffer_d.size());

  std::vector<uint8_t> buffer_a(8);
  a.to_bytes(buffer_a.data());
  std::vector<uint8_t> buffer_b(8);
  b.to_bytes(buffer_b.data());
  REQUIRE(buffer_a == buffer_c);
  REQUIRE(buffer_b == buffer_d);
}
TEST_CASE("NTL to custom conversion", "[field]") {
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

  field::GF2E ab = a * b;
  GF2E a_ntl = utils::custom_to_ntl(a);
  GF2E ab_ntl = utils::custom_to_ntl(ab);
  GF2E b_ntl = ab_ntl / a_ntl;
  field::GF2E b2 = utils::ntl_to_custom(b_ntl);
  REQUIRE(b == b2);
}