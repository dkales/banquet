#pragma once

#include <cstdint>
#include <cstdlib>

class GF2_32 {
  uint64_t data;

public:
  GF2_32() : data(0){};
  GF2_32(uint64_t data) : data(data) {}

  void set_coeff(size_t idx) { data |= (1ULL << idx); }
  GF2_32 operator+(const GF2_32 &other) const;
  GF2_32 operator*(const GF2_32 &other) const;
  bool operator==(const GF2_32 &other) const;

  void to_bytes(uint8_t *out);
  void from_bytes(uint8_t *in);
};