#pragma once

#include "banquet_instances.h"
#include <cstdint>
#include <cstdlib>
#include <functional>
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
  static std::function<uint64_t(__m128i)> reduce;
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

  void to_bytes(uint8_t *out) const;
  std::vector<uint8_t> to_bytes() const;
  void from_bytes(uint8_t *in);
  static void init_extension_field(const banquet_instance_t &instance);

  friend GF2E(::dot_product)(const std::vector<field::GF2E> &lhs,
                             const std::vector<field::GF2E> &rhs);
};

const GF2E &lift_uint8_t(uint8_t value);

std::vector<GF2E> get_first_n_field_elements(size_t n);
std::vector<std::vector<GF2E>>
precompute_lagrange_polynomials(const std::vector<GF2E> &x_values);
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