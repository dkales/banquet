#include "utils.h"

#include <NTL/GF2EX.h>
#include <NTL/GF2X.h>
#include <stdexcept>

using namespace NTL;

namespace utils {

static GF2X modulus;
static std::array<GF2E, 256> lifting_lut;

static void init_lifting_lut(const GF2E &generator) {
  clear(lifting_lut[0]); // lut(0) = 0
  set(lifting_lut[1]);   // lut(1) = 1

  GF2E pow = generator;
  for (size_t bit = 1; bit < 8; bit++) {
    size_t start = (1ULL << bit);
    // copy last half of LUT and add current generator power
    for (size_t idx = 0; idx < start; idx++) {
      lifting_lut[start + idx] = lifting_lut[idx] + pow;
    }
    pow = pow * generator;
  }
}

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
    init_lifting_lut(conv<GF2E>(gen));
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
    init_lifting_lut(conv<GF2E>(gen));
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
    init_lifting_lut(conv<GF2E>(gen));
  } break;
  default:
    throw std::runtime_error(
        "modulus for that specific lambda not implemented.");
  }
}

const GF2E &lift_uint8_t(uint8_t value) { return lifting_lut[value]; }

GF2E GF2E_from_bytes(const std::vector<uint8_t> &value) {
  // assumes value is already smaller than current modulus
  GF2X inner = GF2XFromBytes(value.data(), value.size());
  // GF2E result(INIT_NO_ALLOC);
  // result.LoopHole() = inner;
  // return result;
  return conv<GF2E>(inner);
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

  GF2EX full_poly = BuildFromRoots(x_values);
  GF2EX lagrange_poly;
  GF2EX missing_term;
  SetX(missing_term);
  for (size_t k = 0; k < m; k++) {
    SetCoeff(missing_term, 0, -x_values[k]);
    lagrange_poly = full_poly / missing_term;
    lagrange_poly = lagrange_poly / eval(lagrange_poly, x_values[k]);
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
