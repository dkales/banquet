#include "field.h"

#include <array>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

extern "C" {
#include "portable_endian.h"
}

namespace {
std::array<field::GF2E, 256> lifting_lut;

void init_lifting_lut(const field::GF2E &generator) {
  lifting_lut[0] = field::GF2E(0); // lut(0) = 0
  lifting_lut[1] = field::GF2E(1); // lut(1) = 1

  field::GF2E pow = generator;
  for (size_t bit = 1; bit < 8; bit++) {
    size_t start = (1ULL << bit);
    // copy last half of LUT and add current generator power
    for (size_t idx = 0; idx < start; idx++) {
      lifting_lut[start + idx] = lifting_lut[idx] + pow;
    }
    pow = pow * generator;
  }
}

inline __m128i clmul(uint64_t a, uint64_t b) {
  return _mm_clmulepi64_si128(_mm_set_epi64x(0, a), _mm_set_epi64x(0, b), 0);
}

uint64_t reduce_GF2_16(__m128i in) {
  // modulus = x^16 + x^5 + x^3 + x + 1
  constexpr uint64_t lower_mask = 0xFFFFULL;
  uint64_t R_lower = _mm_cvtsi128_si64(in);
  uint64_t R_upper = R_lower >> 16;

  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 11) ^ (T >> 13) ^ (T >> 15);
  R_lower = R_lower ^ (R_upper << 5) ^ (R_upper << 3) ^ (R_upper << 1) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

// actually a bit slowerthan naive version below
__attribute__((unused)) uint64_t reduce_GF2_32_barret(__m128i in) {
  // modulus = x^32 + x^7 + x^3 + x^2 + 1
  constexpr uint64_t P =
      (1ULL << 32) | (1ULL << 7) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
  constexpr uint64_t mu = P;
  uint64_t R = _mm_cvtsi128_si64(in);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R >> 32, mu));
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1 >> 32, P));
  return 0xFFFFFFFFULL & (R ^ T2);
}
uint64_t reduce_GF2_32(__m128i in) {
  // modulus = x^32 + x^7 + x^3 + x^2 + 1
  constexpr uint64_t lower_mask = 0xFFFFFFFFULL;
  uint64_t R_lower = _mm_cvtsi128_si64(in);
  uint64_t R_upper = R_lower >> 32;

  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 25) ^ (T >> 29) ^ (T >> 30);
  R_lower = R_lower ^ (R_upper << 7) ^ (R_upper << 3) ^ (R_upper << 2) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_40(__m128i in) {
  // modulus = x^40 + x^5 + x^4 + x^3 + 1
  constexpr uint64_t upper_mask = 0xFFFFULL;
  constexpr uint64_t lower_mask = 0xFFFFFFFFFFULL;
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper =
      ((_mm_extract_epi64(in, 1) & upper_mask) << 24) | (R_lower >> 40);

  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 35) ^ (T >> 36) ^ (T >> 37);
  R_lower = R_lower ^ (R_upper << 5) ^ (R_upper << 4) ^ (R_upper << 3) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t reduce_GF2_48(__m128i in) {
  // modulus = x^48 + x^5 + x^3 + x^2 + 1
  constexpr uint64_t upper_mask = 0xFFFFFFFFULL;
  constexpr uint64_t lower_mask = 0xFFFFFFFFFFFFULL;
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper =
      ((_mm_extract_epi64(in, 1) & upper_mask) << 16) | (R_lower >> 48);
  uint64_t T = R_upper;
  R_upper = R_upper ^ (T >> 43) ^ (T >> 45) ^ (T >> 46);
  R_lower = R_lower ^ (R_upper << 5) ^ (R_upper << 3) ^ (R_upper << 2) ^
            (R_upper << 0);
  return lower_mask & R_lower;
}

uint64_t GF2_euclidean_div_quotient(uint64_t a, uint64_t b) {
  uint64_t quotient = 0;
  int diff = __builtin_clzl(b) - __builtin_clzl(a);
  while (diff >= 0 && a != 0) {
    quotient |= (1ULL << diff);
    a ^= (b << diff);
    diff = __builtin_clzl(b) - __builtin_clzl(a);
  }
  return quotient;
}

uint64_t mod_inverse(uint64_t a, uint64_t mod) {
  uint64_t t = 0;
  uint64_t new_t = 1;
  uint64_t r = mod;
  uint64_t new_r = a;
  uint64_t tmp;

  while (new_r != 0) {
    uint64_t quotient = GF2_euclidean_div_quotient(r, new_r);
    tmp = r;
    r = new_r;
    new_r = tmp ^ _mm_extract_epi64(clmul(quotient, new_r), 0);
    tmp = t;
    t = new_t;
    new_t = tmp ^ _mm_extract_epi64(clmul(quotient, new_t), 0);
  }

  return t;
}

} // namespace

namespace field {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
std::function<uint64_t(__m128i)> GF2E::reduce = nullptr;
#pragma GCC diagnostic pop
size_t GF2E::byte_size = 0;
uint64_t GF2E::modulus = 0;

GF2E GF2E::operator+(const GF2E &other) const {
  return GF2E(this->data ^ other.data);
}
GF2E &GF2E::operator+=(const GF2E &other) {
  this->data ^= other.data;
  return *this;
}
GF2E GF2E::operator-(const GF2E &other) const {
  return GF2E(this->data ^ other.data);
}
GF2E &GF2E::operator-=(const GF2E &other) {
  this->data ^= other.data;
  return *this;
}
GF2E GF2E::operator*(const GF2E &other) const {
  return GF2E(reduce(clmul(this->data, other.data)));
}
GF2E &GF2E::operator*=(const GF2E &other) {
  this->data = reduce(clmul(this->data, other.data));
  return *this;
}
bool GF2E::operator==(const GF2E &other) const {
  return this->data == other.data;
}
bool GF2E::operator!=(const GF2E &other) const {
  return this->data != other.data;
}

GF2E GF2E::inverse() const { return GF2E(mod_inverse(this->data, modulus)); }

void GF2E::to_bytes(uint8_t *out) const {
  uint64_t be_data = htole64(data);
  memcpy(out, (uint8_t *)(&be_data), byte_size);
}
std::vector<uint8_t> GF2E::to_bytes() const {
  std::vector<uint8_t> buffer(byte_size);
  this->to_bytes(buffer.data());
  return buffer;
}

void GF2E::from_bytes(uint8_t *in) {
  data = 0;
  memcpy((uint8_t *)(&data), in, byte_size);
  data = le64toh(data);
}

void GF2E::init_extension_field(const banquet_instance_t &instance) {
  switch (instance.lambda) {
  case 2: {
    // modulus = x^16 + x^5 + x^3 + x + 1
    modulus =
        (1ULL << 16) | (1ULL << 5) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
    reduce = reduce_GF2_16;
    byte_size = 2;
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^16
    //   Defn: x |--> y^14 + y^8 + y^7 + y^6 + y^3 + y^2 + 1
    GF2E gen;
    gen.set_coeff(14);
    gen.set_coeff(8);
    gen.set_coeff(7);
    gen.set_coeff(6);
    gen.set_coeff(3);
    gen.set_coeff(2);
    gen.set_coeff(0);

    init_lifting_lut(gen);
  } break;
  case 4: {
    // modulus = x^32 + x^7 + x^3 + x^2 + 1
    modulus =
        (1ULL << 32) | (1ULL << 7) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
    reduce = reduce_GF2_32;
    byte_size = 4;
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^32
    //   Defn: x |--> y^30 + y^23 + y^21 + y^18 + y^14 + y^13 + y^11 + y^9 + y^7
    //   + y^6 + y^5 + y^4 + y^3 + y
    GF2E gen;
    gen.set_coeff(30);
    gen.set_coeff(23);
    gen.set_coeff(21);
    gen.set_coeff(18);
    gen.set_coeff(14);
    gen.set_coeff(13);
    gen.set_coeff(11);
    gen.set_coeff(9);
    gen.set_coeff(7);
    gen.set_coeff(6);
    gen.set_coeff(5);
    gen.set_coeff(4);
    gen.set_coeff(3);
    gen.set_coeff(1);

    init_lifting_lut(gen);
  } break;
  case 5: {
    // modulus = x^40 + x^5 + x^4 + x^3 + 1
    modulus =
        (1ULL << 40) | (1ULL << 5) | (1ULL << 4) | (1ULL << 3) | (1ULL << 0);
    reduce = reduce_GF2_40;
    byte_size = 5;
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^40
    //   Defn: x |--> y^31 + y^30 + y^27 + y^25 + y^22 + y^21 + y^20 + y^18 +
    //   y^15 + y^9 + y^6 + y^4 + y^2
    GF2E gen;
    gen.set_coeff(31);
    gen.set_coeff(30);
    gen.set_coeff(27);
    gen.set_coeff(25);
    gen.set_coeff(22);
    gen.set_coeff(21);
    gen.set_coeff(20);
    gen.set_coeff(18);
    gen.set_coeff(15);
    gen.set_coeff(9);
    gen.set_coeff(6);
    gen.set_coeff(4);
    gen.set_coeff(2);

    init_lifting_lut(gen);
  } break;
  case 6: {
    // modulus = x^48 + x^5 + x^3 + x^2 + 1
    modulus =
        (1ULL << 48) | (1ULL << 5) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
    reduce = reduce_GF2_48;
    byte_size = 6;
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^48
    //   Defn: x |--> y^45 + y^43 + y^40 + y^37 + y^36 + y^35 + y^34 + y^33 +
    //   y^31 + y^30 + y^29 + y^28 + y^24 + y^21 + y^20 + y^19 + y^16 + y^14 +
    //   y^13 + y^11 + y^10 + y^7 + y^3 + y^2
    GF2E gen;
    gen.set_coeff(45);
    gen.set_coeff(43);
    gen.set_coeff(40);
    gen.set_coeff(37);
    gen.set_coeff(36);
    gen.set_coeff(35);
    gen.set_coeff(34);
    gen.set_coeff(33);
    gen.set_coeff(31);
    gen.set_coeff(30);
    gen.set_coeff(29);
    gen.set_coeff(28);
    gen.set_coeff(24);
    gen.set_coeff(21);
    gen.set_coeff(20);
    gen.set_coeff(19);
    gen.set_coeff(16);
    gen.set_coeff(14);
    gen.set_coeff(13);
    gen.set_coeff(11);
    gen.set_coeff(10);
    gen.set_coeff(7);
    gen.set_coeff(3);
    gen.set_coeff(2);

    init_lifting_lut(gen);
  } break;
  default:
    throw std::runtime_error(
        "modulus for that specific lambda not implemented.");
  }
}

const GF2E &lift_uint8_t(uint8_t value) { return lifting_lut[value]; }

std::vector<GF2E> get_first_n_field_elements(size_t n) {
  std::vector<GF2E> result;
  result.reserve(n);
  GF2E x(2);
  GF2E gen = x;
  for (size_t i = 0; i < n; i++) {
    result.push_back(gen);
    gen = gen * x;
  }
  return result;
}

std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials(const std::vector<GF2E> &x_values) {
  size_t m = x_values.size();
  std::vector<std::vector<GF2E>> precomputed_lagrange_polynomials;
  precomputed_lagrange_polynomials.reserve(m);

  std::vector<GF2E> x_except_k;
  GF2E denominator;
  for (size_t k = 0; k < m; k++) {
    denominator = GF2E(1);
    x_except_k.clear();
    x_except_k.reserve(m - 1);
    for (size_t j = 0; j < m; j++) {
      if (k != j) {
        denominator *= x_values[k] - x_values[j];
        x_except_k.push_back(x_values[j]);
      }
    }
    std::vector<GF2E> numerator = build_from_roots(x_except_k);

    numerator = numerator * denominator.inverse();
    precomputed_lagrange_polynomials.push_back(numerator);
  }

  return precomputed_lagrange_polynomials;
}

std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<std::vector<GF2E>> &precomputed_lagrange_polynomials,
    const std::vector<GF2E> &y_values) {
  if (precomputed_lagrange_polynomials.size() != y_values.size() ||
      y_values.empty())
    throw std::runtime_error("invalid sizes for interpolation");

  std::vector<GF2E> res(precomputed_lagrange_polynomials[0].size());
  size_t m = y_values.size();
  for (size_t k = 0; k < m; k++) {
    res += precomputed_lagrange_polynomials[k] * y_values[k];
  }
  return res;
}

std::vector<GF2E> build_from_roots(const std::vector<GF2E> &roots) {
  size_t len = roots.size();

  std::vector<GF2E> poly(roots);
  poly.push_back(GF2E(0));

  GF2E tmp;
  for (size_t k = 1; k < len; k++) {
    tmp = poly[k];
    poly[k] = tmp + poly[k - 1];
    for (size_t i = k - 1; i >= 1; i--) {
      poly[i] = poly[i] * tmp + poly[i - 1];
    }
    poly[0] *= tmp;
  }
  poly[len] = GF2E(1);
  return poly;
}
// horner eval
GF2E eval(const std::vector<GF2E> &poly, const GF2E &point) {
  GF2E acc;
  long i;

  for (i = poly.size() - 1; i >= 0; i--) {
    acc *= point;
    acc += poly[i];
  }

  return acc;
}

} // namespace field

std::vector<field::GF2E> operator+(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs) {
  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");

  std::vector<field::GF2E> result(lhs);
  for (size_t i = 0; i < lhs.size(); i++)
    result[i] += rhs[i];

  return result;
}

std::vector<field::GF2E> &operator+=(std::vector<field::GF2E> &lhs,
                                     const std::vector<field::GF2E> &rhs) {
  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");

  for (size_t i = 0; i < lhs.size(); i++)
    lhs[i] += rhs[i];

  return lhs;
}

// somewhat optimized inner product, only do one lazy reduction
field::GF2E dot_product(const std::vector<field::GF2E> &lhs,
                        const std::vector<field::GF2E> &rhs) {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");

  // field::GF2E result;
  // for (size_t i = 0; i < lhs.size(); i++)
  // result += lhs[i] * rhs[i];
  __m128i accum = _mm_setzero_si128();
  for (size_t i = 0; i < lhs.size(); i++)
    accum = _mm_xor_si128(accum, clmul(lhs[i].data, rhs[i].data));

  field::GF2E result(field::GF2E::reduce(accum));
  return result;
}

std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const field::GF2E &rhs) {
  std::vector<field::GF2E> result(lhs);
  for (size_t i = 0; i < lhs.size(); i++)
    result[i] *= rhs;

  return result;
}

std::vector<field::GF2E> operator*(const field::GF2E &lhs,
                                   const std::vector<field::GF2E> &rhs) {
  return rhs * lhs;
}

// naive polynomial multiplication
std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs) {

  std::vector<field::GF2E> result(lhs.size() + rhs.size() - 1);
  for (size_t i = 0; i < lhs.size(); i++)
    for (size_t j = 0; j < rhs.size(); j++)
      result[i + j] += lhs[i] * rhs[j];

  return result;
}