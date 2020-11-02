#pragma once

#include "banquet.h"
#include <NTL/GF2E.h>
#include <NTL/vec_GF2E.h>
using namespace NTL;

namespace utils {
void init_extension_field(const banquet_instance_t &instance);
GF2E lift_uint8_t(uint8_t value);
GF2E GF2E_from_bytes(const std::vector<uint8_t> &value);

vec_GF2E get_first_n_field_elements(size_t n);
} // namespace utils