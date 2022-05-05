#include "field.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
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

uint64_t reduce_GF2_16_barret(__m128i in) {
  // modulus = x^16 + x^5 + x^3 + x + 1
  constexpr uint64_t lower_mask = 0xFFFFULL;
  constexpr uint64_t P =
      (1ULL << 16) | (1ULL << 5) | (1ULL << 3) | (1ULL << 1) | (1ULL << 0);
  uint64_t R = _mm_cvtsi128_si64(in);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R >> 16, P));
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1 >> 16, P));
  return lower_mask & (R ^ T2);
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
uint64_t reduce_GF2_16_clmul(const __m128i in) {
  // modulus = x^16 + x^5 + x^3 + x + 1
  __m128i p = _mm_set_epi64x(0x0, 0x2B);
  __m128i mask = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFF0000);
  __m128i t;

  __m128i hi = _mm_srli_si128(in, 2); // extracting the in_hi
  __m128i low =
      _mm_xor_si128(_mm_or_si128(in, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // in_hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 4 + 16 -> Length after xor

  hi = _mm_srli_si128(t, 2);                        // extracting the t_hi
  low = _mm_xor_si128(_mm_or_si128(t, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // t_hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 16 -> Length after xor

  return _mm_extract_epi64(t, 0);
}

uint64_t reduce_GF2_32_barret(__m128i in) {
  // modulus = x^32 + x^7 + x^3 + x^2 + 1
  constexpr uint64_t lower_mask = 0xFFFFFFFFULL;
  constexpr uint64_t P =
      (1ULL << 32) | (1ULL << 7) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
  uint64_t R = _mm_cvtsi128_si64(in);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R >> 32, P));
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1 >> 32, P));
  return lower_mask & (R ^ T2);
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
uint64_t reduce_GF2_32_clmul(const __m128i in) {
  // modulus = x^32 + x^7 + x^3 + x^2 + 1
  __m128i p = _mm_set_epi64x(0x0, 0x8d);
  __m128i mask = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFF00000000);
  __m128i t;

  __m128i hi = _mm_srli_si128(in, 4); // extracting the in_hi
  __m128i low =
      _mm_xor_si128(_mm_or_si128(in, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // in_hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 4 + 32 -> Length after xor

  hi = _mm_srli_si128(t, 4);                        // extracting the t_hi
  low = _mm_xor_si128(_mm_or_si128(t, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // t_hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 16 -> Length after xor

  return _mm_extract_epi64(t, 0);
}

uint64_t reduce_GF2_40_barret(__m128i in) {
  // modulus = x^40 + x^5 + x^4 + x^3 + 1
  constexpr uint64_t upper_mask = 0xFFFFULL;
  constexpr uint64_t lower_mask = 0xFFFFFFFFFFULL;
  constexpr uint64_t P =
      (1ULL << 40) | (1ULL << 5) | (1ULL << 4) | (1ULL << 3) | (1ULL << 0);
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper =
      ((_mm_extract_epi64(in, 1) & upper_mask) << 24) | (R_lower >> 40);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R_upper, P));
  T1 = ((_mm_extract_epi64(in, 1) & upper_mask) << 24) | (T1 >> 40);
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1, P));
  return lower_mask & (R_lower ^ T2);
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
uint64_t reduce_GF2_40_clmul(const __m128i in) {
  // modulus = x^40 + x^5 + x^4 + x^3 + 1
  __m128i p = _mm_set_epi64x(0x0, 0x39);
  __m128i mask = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFFFF0000000000);
  __m128i t;

  __m128i hi = _mm_srli_si128(in, 5); // extracting the hi
  __m128i low =
      _mm_xor_si128(_mm_or_si128(in, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 4 + 40 -> Length after xor

  hi = _mm_srli_si128(t, 5);                        // extracting the hi
  low = _mm_xor_si128(_mm_or_si128(t, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 40 -> Length after xor

  return _mm_extract_epi64(t, 0);
}

uint64_t reduce_GF2_48_barret(__m128i in) {
  // modulus = x^48 + x^5 + x^3 + x^2 + 1
  constexpr uint64_t upper_mask = 0xFFFFFFFFULL;
  constexpr uint64_t lower_mask = 0xFFFFFFFFFFFFULL;
  constexpr uint64_t P =
      (1ULL << 48) | (1ULL << 5) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
  uint64_t R_lower = _mm_extract_epi64(in, 0);
  uint64_t R_upper =
      ((_mm_extract_epi64(in, 1) & upper_mask) << 16) | (R_lower >> 48);
  uint64_t T1 = _mm_cvtsi128_si64(clmul(R_upper, P));
  T1 = ((_mm_extract_epi64(in, 1) & upper_mask) << 16) | (T1 >> 48);
  uint64_t T2 = _mm_cvtsi128_si64(clmul(T1, P));
  return lower_mask & (R_lower ^ T2);
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
uint64_t reduce_GF2_48_clmul(const __m128i in) {
  // modulus = x^48 + x^5 + x^3 + x^2 + 1
  __m128i p = _mm_set_epi64x(0x0, 0x2d);
  __m128i mask = _mm_set_epi64x(0xFFFFFFFFFFFFFFFF, 0xFFFF000000000000);
  __m128i t;

  __m128i hi = _mm_srli_si128(in, 6); // extracting the hi
  __m128i low =
      _mm_xor_si128(_mm_or_si128(in, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 4 + 48 -> Length after xor

  hi = _mm_srli_si128(t, 6);                        // extracting the hi
  low = _mm_xor_si128(_mm_or_si128(t, mask), mask); // extracting the low

  t = _mm_clmulepi64_si128(hi, p, 0x00); // hi_low(0x00) * p
  t = _mm_xor_si128(t, low);             // 48 -> Length after xor

  return _mm_extract_epi64(t, 0);
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

uint64_t (*GF2E::reduce_naive)(__m128i);
uint64_t (*GF2E::reduce_barret)(__m128i);
uint64_t (*GF2E::reduce_clmul)(__m128i);

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
  return GF2E(reduce_clmul(clmul(this->data, other.data)));
}
GF2E &GF2E::operator*=(const GF2E &other) {
  this->data = reduce_clmul(clmul(this->data, other.data));
  return *this;
}
bool GF2E::operator==(const GF2E &other) const {
  return this->data == other.data;
}
bool GF2E::operator!=(const GF2E &other) const {
  return this->data != other.data;
}
std::ostream &operator<<(std::ostream &os, const GF2E &ele) {
  os << "0x" << std::setfill('0') << std::hex << std::setw(8) << ele.get_data();
  return os;
}

GF2E GF2E::inverse() const { return GF2E(mod_inverse(this->data, modulus)); }

// This is not faster than calling the operand because we are squaring numbers
// in GF2_X where X is less than 64
GF2E GF2E::sqr() const {

  __m128i tmp[2];
  __m128i res;

  __m128i sqrT = _mm_set_epi64x(0x5554515045444140, 0x1514111005040100);
  __m128i mask = _mm_set_epi64x(0, 0x0F0F0F0F0F0F0F0F);

  __m128i a = _mm_set_epi64x(0, this->data);

  tmp[0] = _mm_and_si128(a, mask);

  tmp[1] = _mm_srli_epi64(a, 4);
  tmp[1] = _mm_and_si128(tmp[1], mask);

  tmp[0] = _mm_shuffle_epi8(sqrT, tmp[0]);
  tmp[1] = _mm_shuffle_epi8(sqrT, tmp[1]);

  res = _mm_unpacklo_epi8(tmp[0], tmp[1]);

  return GF2E(reduce_clmul(res));
}

// This is actually slow than the other inverse inplementation but this has
// constant time, so gets better with larger root size and field size
GF2E GF2E::inverse_const_time() const {

  switch (this->byte_size) {
  case 2: {
    constexpr uint64_t u[6] = {1, 2, 3, 6, 12, 15};
    constexpr uint64_t u_len = sizeof(u) / sizeof(u[0]);
    // q = u[i] - u[i - 1] should give us the corresponding values
    // (1, 1, 3, 6, 3), which will have corresponding indexes
    constexpr uint64_t q_index[u_len - 1] = {0, 0, 2, 3, 2};
    GF2E b[u_len];
    b[0] = this->data;
    for (size_t i = 1; i < u_len; ++i) {
      GF2E b_p = b[i - 1];
      GF2E b_q = b[q_index[i - 1]];
      for (uint64_t m = u[q_index[i - 1]]; m; --m) {
        b_p = b_p.sqr();
      }
      b[i] = reduce_clmul(clmul(b_p.data, b_q.data));
    }
    return b[u_len - 1].sqr();
    break;
  }
  case 4: {
    constexpr uint64_t u[8] = {1, 2, 3, 5, 7, 14, 28, 31};
    constexpr uint64_t u_len = sizeof(u) / sizeof(u[0]);
    // q = u[i] - u[i - 1] should give us the corresponding values
    // (1, 1, 2, 2, 7, 14, 3), which will have corresponding indexes
    constexpr uint64_t q_index[u_len - 1] = {0, 0, 1, 1, 4, 5, 2};
    GF2E b[u_len];
    b[0] = this->data;
    for (size_t i = 1; i < u_len; ++i) {
      GF2E b_p = b[i - 1];
      GF2E b_q = b[q_index[i - 1]];
      for (uint64_t m = u[q_index[i - 1]]; m; --m) {
        b_p = b_p.sqr();
      }
      b[i] = reduce_clmul(clmul(b_p.data, b_q.data));
    }
    return b[u_len - 1].sqr();
    break;
  }
  case 5: {
    constexpr uint64_t u[8] = {1, 2, 3, 6, 12, 24, 27, 39};
    constexpr uint64_t u_len = sizeof(u) / sizeof(u[0]);
    // q = u[i] - u[i - 1] should give us the corresponding values
    // (1, 1, 3, 6, 12, 3, 12), which will have corresponding indexes
    constexpr uint64_t q_index[u_len - 1] = {0, 0, 2, 3, 4, 2, 4};
    GF2E b[u_len];
    b[0] = this->data;
    for (size_t i = 1; i < u_len; ++i) {
      GF2E b_p = b[i - 1];
      GF2E b_q = b[q_index[i - 1]];
      for (uint64_t m = u[q_index[i - 1]]; m; --m) {
        b_p = b_p.sqr();
      }
      b[i] = reduce_clmul(clmul(b_p.data, b_q.data));
    }
    return b[u_len - 1].sqr();
    break;
  }
  case 6: {
    constexpr uint64_t u[9] = {1, 2, 3, 5, 10, 20, 23, 46, 47};
    constexpr uint64_t u_len = sizeof(u) / sizeof(u[0]);
    // q = u[i] - u[i - 1] should give us the corresponding values
    // (1, 1, 2, 5, 10, 3, 23, 1), which will have corresponding indexes
    constexpr uint64_t q_index[u_len - 1] = {0, 0, 1, 3, 4, 2, 6, 0};
    GF2E b[u_len];
    b[0] = this->data;
    for (size_t i = 1; i < u_len; ++i) {
      GF2E b_p = b[i - 1];
      GF2E b_q = b[q_index[i - 1]];
      for (uint64_t m = u[q_index[i - 1]]; m; --m) {
        b_p = b_p.sqr();
      }
      b[i] = reduce_clmul(clmul(b_p.data, b_q.data));
    }
    return b[u_len - 1].sqr();
    break;
  }
  default:
    return GF2E(this->data);
    break;
  }
}

void GF2E::to_bytes(uint8_t *out) const {
  uint64_t be_data = htole64(data);
  memcpy(out, (uint8_t *)(&be_data), byte_size);
}

std::vector<uint8_t> GF2E::to_bytes() const {
  std::vector<uint8_t> buffer(byte_size);
  this->to_bytes(buffer.data());
  return buffer;
}

uint64_t GF2E::get_data() const { return this->data; }

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
    reduce_naive = reduce_GF2_16;
    reduce_barret = reduce_GF2_16_barret;
    reduce_clmul = reduce_GF2_16_clmul;
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
    reduce_naive = reduce_GF2_32;
    reduce_barret = reduce_GF2_32_barret;
    reduce_clmul = reduce_GF2_32_clmul;
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
    reduce_naive = reduce_GF2_40;
    reduce_barret = reduce_GF2_40_barret;
    reduce_clmul = reduce_GF2_40_clmul;
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
    reduce_naive = reduce_GF2_48;
    reduce_barret = reduce_GF2_48_barret;
    reduce_clmul = reduce_GF2_48_clmul;
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

// Use to precompute the constants of the denominaotr.inverse()
std::vector<GF2E> precompute_denominator(const std::vector<GF2E> &x_values) {
  // Check if value size is power of 2
  if (ceil(log2(x_values.size())) != floor(log2(x_values.size()))) {
    throw std::runtime_error("invalid sizes for interpolation");
  }
  size_t values_size = x_values.size();
  std::vector<GF2E> precomputed_denominator;
  precomputed_denominator.reserve(values_size);
  GF2E denominator;

  for (size_t k = 0; k < values_size; ++k) {
    denominator = GF2E(1);
    for (size_t i = 0; i < values_size; ++i) {
      if (i != k) {
        denominator *= x_values[k] - x_values[i];
      }
    }
    precomputed_denominator.push_back(denominator.inverse());
  }

  return precomputed_denominator;
}

void set_x_minus_xi_poly_size(
    std::vector<std::vector<GF2E>> &precomputed_x_minus_xi, size_t root_count) {
  size_t split_count = 0;
  for (size_t i = root_count; i > 1; i /= 2) {
    split_count += i;
  }
  precomputed_x_minus_xi.reserve(split_count);
}

// Use to precompute x - xi recurssively
void precompute_x_minus_xi_poly_splits(
    const std::vector<GF2E> &x_values,
    std::vector<std::vector<GF2E>> &precomputed_x_minus_xi) {

  size_t len = x_values.size();
  if (len == 1) {
    return;
  }
  size_t len_half = len / 2;

  // Gets the first half roots
  std::vector<GF2E> x_first_half_roots;
  x_first_half_roots.reserve(len / 2);
  for (size_t i = 0; i < len_half; ++i) {
    x_first_half_roots.push_back(x_values[i]);
  }
  // Generates poly from roots
  std::vector<GF2E> x_first_half_poly = build_from_roots(x_first_half_roots);
  // Save poly
  precomputed_x_minus_xi.push_back(x_first_half_poly);
  // Recurssion with the first half roots
  precompute_x_minus_xi_poly_splits(x_first_half_roots, precomputed_x_minus_xi);

  // Gets the second half roots
  std::vector<GF2E> x_second_half_roots;
  x_second_half_roots.reserve(len / 2);
  for (size_t i = len_half; i < len; ++i) {
    x_second_half_roots.push_back(x_values[i]);
  }
  // Generates poly from roots
  std::vector<GF2E> x_second_half_poly = build_from_roots(x_second_half_roots);
  // Save poly
  precomputed_x_minus_xi.push_back(x_second_half_poly);
  // Recurssion with the second half roots
  precompute_x_minus_xi_poly_splits(x_second_half_roots,
                                    precomputed_x_minus_xi);
}

// Computing the precomputable part of the plain langrange interpolation
// (not-optimized)
std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials_slow(const std::vector<GF2E> &x_values) {
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

// Computing the precomputable part of the plain langrange interpolation
// (optimized) (fast only with large root sizes)
std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials(const std::vector<GF2E> &x_values) {

  size_t m = x_values.size();
  std::vector<std::vector<GF2E>> precomputed_lagrange_polynomials;
  precomputed_lagrange_polynomials.reserve(m);

  std::vector<GF2E> x_minus_xi = build_from_roots(x_values);
  for (size_t k = 0; k < m; ++k) {

    std::vector<GF2E> numerator = x_minus_xi / x_values[k];

    precomputed_lagrange_polynomials.push_back(
        eval(numerator, x_values[k]).inverse() * numerator);
  }

  return precomputed_lagrange_polynomials;
}

// Langrange interpolation with precomputation (slow)
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

// Langrange interpolation with precomputation (fast)
std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<GF2E> &precomputed_denominator,
    const std::vector<GF2E> &y_values, const size_t index) {

  std::vector<GF2E> a_precomputed_denominator;
  a_precomputed_denominator.reserve(1);
  a_precomputed_denominator.push_back(precomputed_denominator[index]);

  return y_values[index] * a_precomputed_denominator;
}

// Langrange interpolation using recurssion (fast)
std::vector<GF2E> interpolate_with_recurrsion(
    const std::vector<GF2E> &y_values,
    const std::vector<GF2E> &precomputed_denominator,
    const std::vector<std::vector<GF2E>> &precomputed_x_minus_xi,
    const size_t x_start_index, const size_t x_length,
    const size_t x_minus_xi_first_index, const size_t x_minus_xi_length) {

  // For indexing x_values
  const size_t x_len_half = x_length / 2;
  const size_t x_end_index = x_start_index + x_length - 1;
  const size_t x_second_half_start_index = x_start_index + x_len_half;

  // For indexing x_minus_xi values
  const size_t x_minus_xi_half_length = x_minus_xi_length / 2;
  const size_t x_minus_xi_second_index =
      x_minus_xi_half_length + x_minus_xi_first_index;

  // The recurssion part !!
  if (x_length != 2) {

    return (precomputed_x_minus_xi[x_minus_xi_second_index] *
            interpolate_with_recurrsion(y_values, precomputed_denominator,
                                        precomputed_x_minus_xi, x_start_index,
                                        x_len_half, x_minus_xi_first_index + 1,
                                        x_minus_xi_half_length - 1)) +
           (precomputed_x_minus_xi[x_minus_xi_first_index] *
            interpolate_with_recurrsion(
                y_values, precomputed_denominator, precomputed_x_minus_xi,
                x_second_half_start_index, x_len_half,
                x_minus_xi_second_index + 1, x_minus_xi_half_length - 1));
  }

  return (precomputed_x_minus_xi[x_minus_xi_second_index] *
          interpolate_with_precomputation(precomputed_denominator, y_values,
                                          x_start_index)) +
         (precomputed_x_minus_xi[x_minus_xi_first_index] *
          interpolate_with_precomputation(precomputed_denominator, y_values,
                                          x_end_index));
}

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

// normal eval precomputation
std::vector<GF2E> eval_precompute(const GF2E &point, size_t poly_size) {
  std::vector<GF2E> out;
  out.reserve(poly_size);

  GF2E temp = point;

  out.push_back(temp);
  for (size_t i = 1; i < poly_size; ++i) {
    temp *= point;
    out.push_back(temp);
  }
  return out;
}

// normal optmized polynomial evaluation with precomputation optmization
GF2E eval_fast(const std::vector<GF2E> &poly, const std::vector<GF2E> &x_pow_n,
               const size_t lambda) {
  __m128i acc = _mm_set_epi64x(0, poly[0].get_data());
  for (size_t i = 1; i < poly.size(); ++i) {

    acc = acc ^ clmul(poly[i].get_data(), x_pow_n[i - 1].get_data());
  }

  switch (lambda) {
  case 2:
    return GF2E(reduce_GF2_16_clmul(acc));
  case 4:
    return GF2E(reduce_GF2_32_clmul(acc));
  case 5:
    return GF2E(reduce_GF2_40_clmul(acc));
  case 6:
    return GF2E(reduce_GF2_48_clmul(acc));
  default:
    return GF2E(reduce_GF2_32_clmul(acc));
  }
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
    throw std::runtime_error("mul vectors of different sizes");

  __m128i accum = _mm_setzero_si128();
  for (size_t i = 0; i < lhs.size(); i++) {
    accum = _mm_xor_si128(accum, clmul(lhs[i].data, rhs[i].data));
  }
  field::GF2E result(field::GF2E::reduce_clmul(accum));

  return result;
}

// Multiplies polynomial of arbitarty degree
std::vector<field::GF2E>
mul_karatsuba_arbideg(const std::vector<field::GF2E> &lhs,
                      const std::vector<field::GF2E> &rhs) {

  if (lhs.size() != rhs.size())
    throw std::runtime_error("karatsuba mul vectors of different sizes");

  field::GF2E d[lhs.size()];
  std::vector<field::GF2E> c(lhs.size() + rhs.size() - 1);

  // For i == 0
  d[0] = lhs[0] * rhs[0];

  // For i == 1
  d[1] = lhs[1] * rhs[1];
  c[1] = ((lhs[0] + lhs[1]) * (rhs[0] + rhs[1])) - (d[1] + d[0]);

  // For i == 2..poly_length
  for (size_t i = 2; i < lhs.size(); ++i) {
    d[i] = lhs[i] * rhs[i];
    field::GF2E sum;
    for (size_t t = i; t > i / 2; --t) {
      sum +=
          ((lhs[i - t] + lhs[t]) * (rhs[i - t] + rhs[t])) - (d[t] + d[i - t]);
    }
    // If i is even
    sum += d[i / 2] * ((i + 1) % 2);
    c[i] = sum;
  }

  // For i == poly_len..poly_len*2-3
  for (size_t i = lhs.size(); i <= lhs.size() + rhs.size() - 3; ++i) {
    field::GF2E sum;
    for (size_t t = lhs.size() - 1; t > i / 2; --t) {
      sum +=
          ((lhs[i - t] + lhs[t]) * (rhs[i - t] + rhs[t])) - (d[t] + d[i - t]);
    }
    // If i is even
    sum += d[i / 2] * ((i + 1) % 2);
    c[i] = sum;
  }

  // Setting the first and the last i
  c[0] = d[0];
  c[lhs.size() + rhs.size() - 2] = d[lhs.size() - 1];

  return c;
}

// Adding dummy values to make the coeff size a power of 2
void mul_karatsuba_fixdeg_precondition_poly(std::vector<field::GF2E> &lhs,
                                            std::vector<field::GF2E> &rhs) {

  // Computes the next power of 2 for a 32 bit number
  // This represents the no. of coeffs
  size_t next_2_pow = lhs.size();
  next_2_pow--;
  next_2_pow |= next_2_pow >> 1;
  next_2_pow |= next_2_pow >> 2;
  next_2_pow |= next_2_pow >> 4;
  next_2_pow |= next_2_pow >> 8;
  next_2_pow |= next_2_pow >> 16;
  next_2_pow++;

  // Putting the dummy terms in the begining to make the polys 2^n - 1 degree ->
  // 2^n coeff
  lhs.resize(next_2_pow, field::GF2E(0));
  rhs.resize(next_2_pow, field::GF2E(0));
}

// Adding dummy values to make the coeff size a power of 2
void mul_karatsuba_fixdeg_normalize_poly(std::vector<field::GF2E> &poly,
                                         size_t old_size) {
  poly.resize((old_size << 1) - 1);
}

// Multiplies polynomial of 2^n - 1 degree -> 2^n coeff
std::vector<field::GF2E>
mul_karatsuba_fixdeg(const std::vector<field::GF2E> &lhs,
                     const std::vector<field::GF2E> &rhs,
                     const size_t start_idx, const size_t end_idx) {

  size_t full_size = ((end_idx - start_idx) + 1);
  size_t half_size = full_size / 2;

  // If a polynomial with degree 0 -> const
  if (full_size == 1) {
    return std::vector<field::GF2E>(1, lhs[start_idx] * rhs[start_idx]);
  }

  std::vector<field::GF2E> d_0 =
      mul_karatsuba_fixdeg(lhs, rhs, start_idx, (start_idx + half_size) - 1);
  std::vector<field::GF2E> d_1 =
      mul_karatsuba_fixdeg(lhs, rhs, (start_idx + half_size), end_idx);

  std::vector<field::GF2E> lhs_l_add_u, rhs_l_add_u;
  lhs_l_add_u.reserve(half_size);
  rhs_l_add_u.reserve(half_size);
  // Getting the lower of lhs and rhs
  for (size_t i = start_idx; i < start_idx + half_size; ++i) {
    lhs_l_add_u.push_back(lhs[i] + lhs[half_size + i]);
    rhs_l_add_u.push_back(rhs[i] + rhs[half_size + i]);
  }
  std::vector<field::GF2E> d_01 =
      mul_karatsuba_fixdeg(lhs_l_add_u, rhs_l_add_u, 0, half_size - 1);

  std::vector<field::GF2E> c(d_1.size() + full_size);
  // D_1*x^n + (D_01 - D_0 - D_1)*x^(n/2) + d_0
  size_t cidx = c.size();
  for (size_t i = d_1.size(); i; --i) {
    c[cidx - 1] += d_1[i - 1];
    cidx--;
  }
  cidx = c.size() * 3 / 4;
  for (size_t i = d_0.size(); i; --i) {
    c[cidx - 1] += (d_01[i - 1] - d_0[i - 1] - d_1[i - 1]);
    cidx--;
  }
  cidx = (c.size() / 2);
  for (size_t i = d_0.size(); i; --i) {
    c[cidx - 1] += d_0[i - 1];
    cidx--;
  }

  return c;
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

  // return mul_karatsuba_arbideg(lhs, rhs);
}

// polynomial division of x^n*c^n + x^n-1*c^n-1 + .... by x - a
// -> Returns Quotient
std::vector<field::GF2E> operator/(const std::vector<field::GF2E> &lhs,
                                   const field::GF2E &rhs) {
  std::vector<field::GF2E> temp(lhs);
  size_t end_index = temp.size() - 1;

  std::vector<field::GF2E> result;
  result.reserve(end_index);

  for (size_t i = end_index; i; --i) {
    temp[i - 1] -= temp[i] * rhs;
  }
  for (size_t i = 1; i <= end_index; ++i) {
    result.push_back(temp[i]);
  }

  return result;
}
