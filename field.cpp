#include "field.h"

#include <cstring>
#include <stdexcept>
extern "C" {
#include "endian_compat.h"
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
  return _mm_clmulepi64_si128(_mm_set1_epi64x(a), _mm_set1_epi64x(b), 0);
}

uint64_t reduce_GF2_32(__m128i in) {
  // modulus = x^32 + x^7 + x^3 + x^2 + 1
  constexpr uint64_t P =
      (1ULL << 32) | (1ULL << 7) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
  constexpr uint64_t mu = P;
  uint64_t R = _mm_extract_epi64(in, 0);
  uint64_t T1 = _mm_extract_epi64(clmul(R >> 32, mu), 0);
  uint64_t T2 = _mm_extract_epi64(clmul(T1 >> 32, P), 0);
  return 0xFFFFFFFFULL & (R ^ T2);
}
uint64_t reduce_GF2_40(__m128i in) {
  // modulus = x^40 + x^5 + x^4 + x^3 + 1
  constexpr uint64_t P =
      (1ULL << 40) | (1ULL << 5) | (1ULL << 4) | (1ULL << 3) | (1ULL << 0);
  constexpr uint64_t mu = P;
  constexpr uint64_t upper_mask = 0xFFFFULL;
  uint64_t R = _mm_extract_epi64(in, 0);
  uint64_t R_red = ((_mm_extract_epi64(in, 1) & upper_mask) << 24) | R >> 40;
  __m128i T1 = clmul(R_red, mu);
  uint64_t T1_red = ((_mm_extract_epi64(T1, 1) & upper_mask) << 24) |
                    _mm_extract_epi64(T1, 0) >> 40;
  uint64_t T2 = _mm_extract_epi64(clmul(T1_red, P), 0);
  return 0xFFFFFFFFFFULL & (R ^ T2);
}

uint64_t reduce_GF2_48(__m128i in) {
  // modulus = x^48 + x^5 + x^3 + x^2 + 1
  constexpr uint64_t P =
      (1ULL << 48) | (1ULL << 5) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
  constexpr uint64_t mu = P;
  constexpr uint64_t upper_mask = 0xFFFFFFFFULL;
  uint64_t R = _mm_extract_epi64(in, 0);
  uint64_t R_red = ((_mm_extract_epi64(in, 1) & upper_mask) << 16) | R >> 48;
  __m128i T1 = clmul(R_red, mu);
  uint64_t T1_red = ((_mm_extract_epi64(T1, 1) & upper_mask) << 16) |
                    _mm_extract_epi64(T1, 0) >> 48;
  uint64_t T2 = _mm_extract_epi64(clmul(T1_red, P), 0);
  return 0xFFFFFFFFFFFFFFULL & (R ^ T2);
}

} // namespace

namespace field {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
std::function<uint64_t(__m128i)> GF2E::reduce = reduce_GF2_32;
#pragma GCC diagnostic pop
size_t GF2E::byte_size = 4;

GF2E GF2E::operator+(const GF2E &other) const {
  return GF2E(this->data ^ other.data);
}
GF2E &GF2E::operator+=(const GF2E &other) {
  this->data ^= other.data;
  return *this;
}
GF2E GF2E::operator*(const GF2E &other) const {
  return GF2E(reduce(clmul(this->data, other.data)));
}
bool GF2E::operator==(const GF2E &other) const {
  return this->data == other.data;
}

void GF2E::to_bytes(uint8_t *out) const {
  uint64_t be_data = htole64(data);
  memcpy(out, (uint8_t *)(&be_data), byte_size);
}
void GF2E::from_bytes(uint8_t *in) {
  data = 0;
  memcpy((uint8_t *)(&data), in, byte_size);
  data = le64toh(data);
}

void GF2E::init_extension_field(const banquet_instance_t &instance) {
  switch (instance.lambda) {
  case 4: {
    // modulus = x^32 + x^7 + x^3 + x^2 + 1
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

  // std::vector<GF2E> full_poly = build_from_roots(x_values);
  std::vector<GF2E> lagrange_poly;
  std::vector<GF2E> missing_term;
  // SetX(missing_term);
  for (size_t k = 0; k < m; k++) {
    // SetCoeff(missing_term, 0, -x_values[k]);
    // lagrange_poly = full_poly / missing_term;
    // lagrange_poly = lagrange_poly / eval(lagrange_poly, x_values[k]);
    precomputed_lagrange_polynomials.push_back(lagrange_poly);
  }

  return precomputed_lagrange_polynomials;
}

std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<std::vector<GF2E>> &precomputed_lagrange_polynomials,
    const std::vector<GF2E> &y_values) {
  if (precomputed_lagrange_polynomials.size() != y_values.size())
    throw std::runtime_error("invalid sizes for interpolation");

  std::vector<GF2E> res;
  size_t m = y_values.size();
  for (size_t k = 0; k < m; k++) {
    // todo vector * const mul
    // res += precomputed_lagrange_polynomials[k] * y_values[k];
  }
  return res;
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
field::GF2E operator*(const std::vector<field::GF2E> &lhs,
                      const std::vector<field::GF2E> &rhs) {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("adding vectors of different sizes");
  __m128i accum = _mm_setzero_si128();
  for (size_t i = 0; i < lhs.size(); i++)
    accum = _mm_xor_si128(accum, clmul(lhs[i].data, rhs[i].data));

  field::GF2E result(field::GF2E::reduce(accum));
  return result;
}