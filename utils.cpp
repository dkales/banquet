#include "utils.h"

#include <NTL/GF2EX.h>
#include <NTL/GF2X.h>
#include <stdexcept>

using namespace NTL;

namespace utils {

static GF2X modulus;
static GF2X generator;

void init_extension_field(const banquet_instance_t &instance) {
  switch (instance.lambda) {
  case 4: {
    // modulus = x^32 + x^7 + x^3 + x^2 + 1
    clear(modulus);
    SetCoeff(modulus, 32);
    SetCoeff(modulus, 7);
    SetCoeff(modulus, 3);
    SetCoeff(modulus, 2);
    SetCoeff(modulus, 0);
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^32
    //   Defn: x |--> y^30 + y^23 + y^21 + y^18 + y^14 + y^13 + y^11 + y^9 + y^7
    //   + y^6 + y^5 + y^4 + y^3 + y
    clear(generator);
    SetCoeff(generator, 30);
    SetCoeff(generator, 23);
    SetCoeff(generator, 21);
    SetCoeff(generator, 18);
    SetCoeff(generator, 14);
    SetCoeff(generator, 13);
    SetCoeff(generator, 11);
    SetCoeff(generator, 9);
    SetCoeff(generator, 7);
    SetCoeff(generator, 6);
    SetCoeff(generator, 5);
    SetCoeff(generator, 4);
    SetCoeff(generator, 3);
    SetCoeff(generator, 1);

    GF2E::init(modulus);
  } break;
  case 5: {
    // modulus = x^40 + x^5 + x^4 + x^3 + 1
    clear(modulus);
    SetCoeff(modulus, 40);
    SetCoeff(modulus, 5);
    SetCoeff(modulus, 4);
    SetCoeff(modulus, 3);
    SetCoeff(modulus, 0);
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^40
    //   Defn: x |--> y^31 + y^30 + y^27 + y^25 + y^22 + y^21 + y^20 + y^18 +
    //   y^15 + y^9 + y^6 + y^4 + y^2
    clear(generator);
    SetCoeff(generator, 31);
    SetCoeff(generator, 30);
    SetCoeff(generator, 27);
    SetCoeff(generator, 25);
    SetCoeff(generator, 22);
    SetCoeff(generator, 21);
    SetCoeff(generator, 20);
    SetCoeff(generator, 18);
    SetCoeff(generator, 15);
    SetCoeff(generator, 9);
    SetCoeff(generator, 6);
    SetCoeff(generator, 4);
    SetCoeff(generator, 2);

    GF2E::init(modulus);
  } break;
  case 6: {
    // modulus = x^48 + x^5 + x^3 + x^2 + 1
    clear(modulus);
    SetCoeff(modulus, 48);
    SetCoeff(modulus, 5);
    SetCoeff(modulus, 3);
    SetCoeff(modulus, 2);
    SetCoeff(modulus, 0);
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^48
    //   Defn: x |--> y^45 + y^43 + y^40 + y^37 + y^36 + y^35 + y^34 + y^33 +
    //   y^31 + y^30 + y^29 + y^28 + y^24 + y^21 + y^20 + y^19 + y^16 + y^14 +
    //   y^13 + y^11 + y^10 + y^7 + y^3 + y^2
    clear(generator);
    SetCoeff(generator, 45);
    SetCoeff(generator, 43);
    SetCoeff(generator, 40);
    SetCoeff(generator, 37);
    SetCoeff(generator, 36);
    SetCoeff(generator, 35);
    SetCoeff(generator, 34);
    SetCoeff(generator, 33);
    SetCoeff(generator, 31);
    SetCoeff(generator, 30);
    SetCoeff(generator, 29);
    SetCoeff(generator, 28);
    SetCoeff(generator, 24);
    SetCoeff(generator, 21);
    SetCoeff(generator, 20);
    SetCoeff(generator, 19);
    SetCoeff(generator, 16);
    SetCoeff(generator, 14);
    SetCoeff(generator, 13);
    SetCoeff(generator, 11);
    SetCoeff(generator, 10);
    SetCoeff(generator, 7);
    SetCoeff(generator, 3);
    SetCoeff(generator, 2);

    GF2E::init(modulus);
  } break;
  default:
    throw std::runtime_error(
        "modulus for that specific lambda not implemented.");
  }
}

GF2E lift_uint8_t(uint8_t value) {

  GF2X result = GF2XFromBytes(&value, sizeof(uint8_t));
  result *= generator;
  return conv<GF2E>(result);
}

GF2E GF2E_from_bytes(const std::vector<uint8_t> &value) {
  GF2X res = GF2XFromBytes(value.data(), value.size());
  return conv<GF2E>(res);
}

vec_GF2E get_first_n_field_elements(size_t n) {
  vec_GF2E result;
  result.SetLength(n);
  GF2X gen;
  SetX(gen);
  for (size_t i = 0; i < n; i++) {
    result[i] = conv<GF2E>(gen);
    gen = MulByX(gen);
  }
  return result;
}
std::vector<GF2EX> precompute_lagrange_polynomials(const vec_GF2E &x_values) {
  size_t m = x_values.length();
  std::vector<GF2EX> precomputed_lagrange_polynomials;
  precomputed_lagrange_polynomials.reserve(m);

  for (size_t k = 0; k < m; k++) {
    vec_GF2E x_without_xk;
    GF2E denom;
    set(denom); // denom = 1
    for (size_t t = 0; t < m; t++) {
      if (t != k) {
        x_without_xk.append(x_values[t]);
        denom *= (x_values[k] - x_values[t]);
      }
    }
    GF2EX lagrange_poly = BuildFromRoots(x_without_xk);
    lagrange_poly = lagrange_poly / denom;
    precomputed_lagrange_polynomials.push_back(lagrange_poly);
  }

  return precomputed_lagrange_polynomials;
}

GF2EX interpolate_with_precomputation(
    const std::vector<GF2EX> &precomputed_lagrange_polynomials,
    const vec_GF2E &y_values) {
  if (precomputed_lagrange_polynomials.size() != y_values.length())
    throw std::runtime_error("invalid sizes for interpolation");

  GF2EX res;
  size_t m = y_values.length();
  for (size_t k = 0; k < m; k++) {
    res += precomputed_lagrange_polynomials[k] * y_values[k];
  }
  return res;
}
} // namespace utils