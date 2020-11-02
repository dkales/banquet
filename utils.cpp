#include "utils.h"

#include <NTL/GF2X.h>
#include <stdexcept>

using namespace NTL;

namespace utils {

static GF2X modulus;
static GF2X generator;

void init_GF2E_modulus(const banquet_instance_t &instance) {
  switch (instance.lambda) {
  case 4: {
    // modulus = X^32 + X^22 + X^2 + X^1 + 1
    clear(modulus);
    SetCoeff(modulus, 32);
    SetCoeff(modulus, 22);
    SetCoeff(modulus, 2);
    SetCoeff(modulus, 1);
    SetCoeff(modulus, 0);
    // Ring morphism:
    //   From: Finite Field in x of size 2^8
    //   To:   Finite Field in y of size 2^32
    //   Defn: x |--> y^27 + y^23 + y^22 + y^21 + y^19 + y^17 + y^13 + y^12 +
    //   y^6
    //   + y^5 + y^2 + 1
    clear(generator);
    SetCoeff(generator, 27);
    SetCoeff(generator, 23);
    SetCoeff(generator, 22);
    SetCoeff(generator, 21);
    SetCoeff(generator, 19);
    SetCoeff(generator, 17);
    SetCoeff(generator, 13);
    SetCoeff(generator, 12);
    SetCoeff(generator, 6);
    SetCoeff(generator, 5);
    SetCoeff(generator, 2);
    SetCoeff(generator, 0);

    GF2E::init(modulus);
  } break;
  default:
    throw std::runtime_error(
        "modulus for that specific lambda not implemented.");
  }
}

GF2E lift_uint8_t(uint8_t value) {

  GF2X result = GF2XFromBytes(&value, sizeof(uint8_t));
  result *= generator;
  return conv<GF2E>(result);
}

GF2E GF2E_from_bytes(const std::vector<uint8_t> &value) {
  GF2X res = GF2XFromBytes(value.data(), value.size());
  return conv<GF2E>(res);
}

vec_GF2E get_first_n_field_elements(size_t n) {
  vec_GF2E result;
  result.SetLength(n);
  GF2X gen;
  SetX(gen);
  for (size_t i = 0; i < n; n++) {
    result[i] = conv<GF2E>(gen);
    gen = MulByX(gen);
  }
  return result;
}
} // namespace utils