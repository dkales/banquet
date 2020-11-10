#include "field.h"

#include <cstring>
extern "C" {
#include "endian_compat.h"
#include <smmintrin.h>
#include <wmmintrin.h>
}

namespace {
inline __m128i clmul(uint64_t a, uint64_t b) {
  return _mm_clmulepi64_si128(_mm_set1_epi64x(a), _mm_set1_epi64x(b), 0);
}

uint64_t reduce_GF2_32(__m128i in) {
  uint64_t P =
      (1ULL << 32) | (1ULL << 7) | (1ULL << 3) | (1ULL << 2) | (1ULL << 0);
  uint64_t mu = P;
  uint64_t R = _mm_extract_epi64(in, 0);
  uint64_t T1 = _mm_extract_epi64(clmul(R >> 32, mu), 0);
  uint64_t T2 = _mm_extract_epi64(clmul(T1 >> 32, P), 0);
  return 0xFFFFFFFFULL & (R ^ T2);
}

} // namespace
GF2_32 GF2_32::operator+(const GF2_32 &other) const {
  return GF2_32(this->data ^ other.data);
}
GF2_32 GF2_32::operator*(const GF2_32 &other) const {
  return GF2_32(reduce_GF2_32(clmul(this->data, other.data)));
}
bool GF2_32::operator==(const GF2_32 &other) const {
  return this->data == other.data;
}

void GF2_32::to_bytes(uint8_t *out) {
  uint64_t be_data = htobe64(data);
  memcpy(out, (uint8_t *)(&be_data), sizeof(be_data));
}
void GF2_32::from_bytes(uint8_t *in) {
  memcpy((uint8_t *)(&data), in, sizeof(data));
  data = be64toh(data);
}