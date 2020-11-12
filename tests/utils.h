#pragma once

#include "../banquet.h"
#include "../field.h"
#include <NTL/GF2E.h>
#include <NTL/vec_GF2E.h>

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

} // namespace utils