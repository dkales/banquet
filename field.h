#pragma once

#include "banquet_instances.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
extern "C" {
#include <smmintrin.h>
#include <wmmintrin.h>
}

namespace field {
class GF2E;
}

field::GF2E dot_product(const std::vector<field::GF2E> &lhs,
                        const std::vector<field::GF2E> &rhs);

namespace field {
class GF2E {

  uint64_t data;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
  static uint64_t (*reduce_naive)(__m128i);
  static uint64_t (*reduce_barret)(__m128i);
  static uint64_t (*reduce_clmul)(__m128i);
#pragma GCC diagnostic pop
  static size_t byte_size;
  static uint64_t modulus;

public:
  GF2E() : data(0){};
  GF2E(uint64_t data) : data(data) {}
  GF2E(const GF2E &other) = default;
  ~GF2E() = default;
  GF2E &operator=(const GF2E &other) = default;

  void clear() { data = 0; }
  void set_coeff(size_t idx) { data |= (1ULL << idx); }
  GF2E operator+(const GF2E &other) const;
  GF2E &operator+=(const GF2E &other);
  GF2E operator-(const GF2E &other) const;
  GF2E &operator-=(const GF2E &other);
  GF2E operator*(const GF2E &other) const;
  GF2E &operator*=(const GF2E &other);
  bool operator==(const GF2E &other) const;
  bool operator!=(const GF2E &other) const;

  GF2E inverse() const;

  GF2E sqr() const;

  GF2E inverse_const_time() const;

  void to_bytes(uint8_t *out) const;
  std::vector<uint8_t> to_bytes() const;
  void from_bytes(uint8_t *in);
  static void init_extension_field(const banquet_instance_t &instance);

  friend GF2E(::dot_product)(const std::vector<field::GF2E> &lhs,
                             const std::vector<field::GF2E> &rhs);

  GF2E(std::string hex_string) {
    // check if hex_string start with 0x or 0X
    if (hex_string.rfind("0x", 0) == 0 || hex_string.rfind("0X", 0) == 0) {
      hex_string = hex_string.substr(2);
    } else {
      throw std::runtime_error("input needs to be a hex number");
    }
    constexpr size_t num_hex_chars = 64 / 4;
    if (hex_string.length() > num_hex_chars)
      throw std::runtime_error("input hex is too large");
    // pad to 128 bit
    hex_string.insert(hex_string.begin(), num_hex_chars - hex_string.length(),
                      '0');

    data = std::stoull(hex_string.substr(0, 64 / 4), nullptr, 16);
  }

  uint64_t get_data() const;
};

std::ostream &operator<<(std::ostream &os, const GF2E &ele);

const GF2E &lift_uint8_t(uint8_t value);

std::vector<GF2E> precompute_denominator(const std::vector<GF2E> &x_values);

void set_x_minus_xi_poly_size(
    std::vector<std::vector<GF2E>> &precomputed_x_minus_xi, size_t root_count);

void precompute_x_minus_xi_poly_splits(
    const std::vector<GF2E> &x_values,
    std::vector<std::vector<GF2E>> &precomputed_x_minus_xi);

std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials_slow(const std::vector<GF2E> &x_values);

std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials(const std::vector<GF2E> &x_values);

std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<std::vector<GF2E>> &precomputed_lagrange_polynomials,
    const std::vector<GF2E> &y_values);

std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<GF2E> &precomputed_denominator,
    const std::vector<GF2E> &y_values, const size_t index);

std::vector<GF2E> interpolate_with_recurrsion(
    const std::vector<GF2E> &y_values,
    const std::vector<GF2E> &precomputed_denominator,
    const std::vector<std::vector<GF2E>> &precomputed_x_minus_xi,
    const size_t x_start_index, const size_t x_length,
    const size_t x_minus_xi_first_index, const size_t x_minus_xi_length);

std::vector<GF2E> get_first_n_field_elements(size_t n);

std::vector<GF2E> interpolate_with_precomputation(
    const std::vector<std::vector<GF2E>> &precomputed_lagrange_polynomials,
    const std::vector<GF2E> &y_values);

std::vector<GF2E> build_from_roots(const std::vector<GF2E> &roots);

GF2E eval(const std::vector<GF2E> &poly, const GF2E &point);

} // namespace field

std::vector<field::GF2E> operator+(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> &operator+=(std::vector<field::GF2E> &self,
                                     const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const field::GF2E &rhs);
std::vector<field::GF2E> operator*(const field::GF2E &lhs,
                                   const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> operator*(const std::vector<field::GF2E> &lhs,
                                   const std::vector<field::GF2E> &rhs);
std::vector<field::GF2E> operator/(const std::vector<field::GF2E> &lhs,
                                   const field::GF2E &rhs);