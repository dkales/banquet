#pragma once

#include "banquet.h"
#include "field.h"
#include <NTL/GF2E.h>
#include <NTL/vec_GF2E.h>
#include <span>

using namespace NTL;

namespace utils {
field::GF2E ntl_to_custom(const GF2E &element);
GF2E custom_to_ntl(const field::GF2E &element);
void init_extension_field(const banquet_instance_t &instance);
const GF2E &lift_uint8_t(uint8_t value);
GF2E GF2E_from_bytes(const std::vector<uint8_t> &value);

vec_GF2E get_first_n_field_elements(size_t n);
std::vector<GF2EX> precompute_lagrange_polynomials(const vec_GF2E &x_values);
GF2EX interpolate_with_precomputation(
    const std::vector<GF2EX> &precomputed_lagrange_polynomials,
    const vec_GF2E &y_values);

class RepByteContainer {
  std::vector<uint8_t> _data;
  size_t _num_repetitions;
  size_t _num_parties;
  size_t _object_size;

public:
  RepByteContainer(size_t num_repetitions, size_t num_parties,
                   size_t object_size)
      : _data(num_repetitions * num_parties * object_size),
        _num_repetitions(num_repetitions), _num_parties(num_parties),
        _object_size(object_size) {}

  inline std::span<uint8_t> get(size_t repetition, size_t party) {
    size_t offset =
        (repetition * _num_parties * _object_size) + (party * _object_size);
    return std::span<uint8_t>(_data.data() + offset, _object_size);
  }

  std::vector<std::span<uint8_t>> get_repetition(size_t repetition) {
    std::vector<std::span<uint8_t>> ret;
    ret.reserve(_num_parties);
    size_t offset = (repetition * _num_parties * _object_size);
    for (size_t i = 0; i < _num_parties; i++)
      ret.emplace_back(_data.data() + offset + i * _object_size, _object_size);
    return ret;
  }
};
} // namespace utils