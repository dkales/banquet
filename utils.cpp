#include "utils.h"

#include <NTL/GF2EX.h>
#include <NTL/GF2X.h>
#include <stdexcept>

using namespace NTL;

namespace utils {

static GF2X modulus;
static std::array<GF2E, 8> generator_powers;

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
    GF2X gen;
    clear(gen);
    SetCoeff(gen, 30);
    SetCoeff(gen, 23);
    SetCoeff(gen, 21);
    SetCoeff(gen, 18);
    SetCoeff(gen, 14);
    SetCoeff(gen, 13);
    SetCoeff(gen, 11);
    SetCoeff(gen, 9);
    SetCoeff(gen, 7);
    SetCoeff(gen, 6);
    SetCoeff(gen, 5);
    SetCoeff(gen, 4);
    SetCoeff(gen, 3);
    SetCoeff(gen, 1);

    GF2E::init(modulus);
    set(generator_powers[0]);
    for (size_t i = 1; i < 8; i++)
      generator_powers[i] = generator_powers[i - 1] * conv<GF2E>(gen);
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
    GF2X gen;
    clear(gen);
    SetCoeff(gen, 31);
    SetCoeff(gen, 30);
    SetCoeff(gen, 27);
    SetCoeff(gen, 25);
    SetCoeff(gen, 22);
    SetCoeff(gen, 21);
    SetCoeff(gen, 20);
    SetCoeff(gen, 18);
    SetCoeff(gen, 15);
    SetCoeff(gen, 9);
    SetCoeff(gen, 6);
    SetCoeff(gen, 4);
    SetCoeff(gen, 2);

    GF2E::init(modulus);
    set(generator_powers[0]);
    for (size_t i = 1; i < 8; i++)
      generator_powers[i] = generator_powers[i - 1] * conv<GF2E>(gen);
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
    GF2X gen;
    clear(gen);
    SetCoeff(gen, 45);
    SetCoeff(gen, 43);
    SetCoeff(gen, 40);
    SetCoeff(gen, 37);
    SetCoeff(gen, 36);
    SetCoeff(gen, 35);
    SetCoeff(gen, 34);
    SetCoeff(gen, 33);
    SetCoeff(gen, 31);
    SetCoeff(gen, 30);
    SetCoeff(gen, 29);
    SetCoeff(gen, 28);
    SetCoeff(gen, 24);
    SetCoeff(gen, 21);
    SetCoeff(gen, 20);
    SetCoeff(gen, 19);
    SetCoeff(gen, 16);
    SetCoeff(gen, 14);
    SetCoeff(gen, 13);
    SetCoeff(gen, 11);
    SetCoeff(gen, 10);
    SetCoeff(gen, 7);
    SetCoeff(gen, 3);
    SetCoeff(gen, 2);

    GF2E::init(modulus);
    set(generator_powers[0]);
    for (size_t i = 1; i < 8; i++)
      generator_powers[i] = generator_powers[i - 1] * conv<GF2E>(gen);
  } break;
  default:
    throw std::runtime_error(
        "modulus for that specific lambda not implemented.");
  }
}

#if 0
std::vector<GF2E> g_precomputed_lifts(256);
GF2E lift_uint8_t(uint8_t value) {
  GF2E result;

    if(g_precomputed_lifts[value] != 0) {
        return g_precomputed_lifts[value];
    }

  for (size_t bit = 0; bit < 8; bit++) {
    GF2 value_bit = conv<GF2>((value >> bit) & 1);
    result += value_bit * generator_powers[bit];
  }

  g_precomputed_lifts[value] = result;
  return result;
}
#else
GF2E lift_uint8_t(uint8_t value) {
  GF2E result;

  for (size_t bit = 0; bit < 8; bit++) {
    GF2 value_bit = conv<GF2>((value >> bit) & 1);
    result += value_bit * generator_powers[bit];
  }

  return result;
}
#endif

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
