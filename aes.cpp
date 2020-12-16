#include "aes.h"
#include <cassert>
extern "C" {
#include <wmmintrin.h> //for intrinsics for AES-NI
}

namespace {

#define ROTL8(x, shift) ((uint8_t)((x) << (shift)) | ((x) >> (8 - (shift))))
constexpr unsigned char AES_SBOX_AFFINE_CONST = 0x63;

unsigned char multiply(unsigned int a, unsigned int b) {
  unsigned char result = 0;

  for (int i = 0; i < 8; ++i) {
    uint8_t mask = -(b & 1);
    result ^= (a & mask);
    uint16_t mask2 = -((uint16_t)(a >> 7) & 1);
    a <<= 1;
    a ^= (0x11b & mask2);
    b >>= 1;
  }
  return result;
}

unsigned char square(unsigned char c) { return multiply(c, c); }

unsigned char xtime(unsigned char c) { return multiply(c, 2); }

unsigned char bytesub(unsigned char c) {
  unsigned char c3 = multiply(square(c), c);
  unsigned char c7 = multiply(square(c3), c);
  unsigned char c63 = multiply(square(square(square(c7))), c7);
  unsigned char c127 = multiply(square(c63), c);
  unsigned char c254 = square(c127);

  unsigned char result = c254 ^ ROTL8(c254, 1) ^ ROTL8(c254, 2) ^
                         ROTL8(c254, 3) ^ ROTL8(c254, 4) ^
                         AES_SBOX_AFFINE_CONST;

  return result;
}
unsigned char inverse(unsigned char c) {
  unsigned char c3 = multiply(square(c), c);
  unsigned char c7 = multiply(square(c3), c);
  unsigned char c63 = multiply(square(square(square(c7))), c7);
  unsigned char c127 = multiply(square(c63), c);
  unsigned char c254 = square(c127);
  return c254;
}
unsigned char bytesub_save(unsigned char c, std::pair<uint8_t, uint8_t> &save) {
  unsigned char c3 = multiply(square(c), c);
  unsigned char c7 = multiply(square(c3), c);
  unsigned char c63 = multiply(square(square(square(c7))), c7);
  unsigned char c127 = multiply(square(c63), c);
  unsigned char c254 = square(c127);
  save.first = c;
  save.second = c254;

  unsigned char result = c254 ^ ROTL8(c254, 1) ^ ROTL8(c254, 2) ^
                         ROTL8(c254, 3) ^ ROTL8(c254, 4) ^
                         AES_SBOX_AFFINE_CONST;
  return result;
}
// does not add the affine constant, this is added to first party manually
inline unsigned char bytesub_restore(unsigned char t) {
  unsigned char result =
      t ^ ROTL8(t, 1) ^ ROTL8(t, 2) ^ ROTL8(t, 3) ^ ROTL8(t, 4);
  return result;
}
// AES-NI functions
#define AES_128_key_exp(k, rcon)                                               \
  aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

static __m128i aes_128_key_expansion(__m128i key, __m128i keygened) {
  keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3, 3, 3, 3));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  return _mm_xor_si128(key, keygened);
}
#define AES_128_key_exp_restore(key, m, s_shares, t_shares, party, sbox_index, \
                                affine, rcon)                                  \
  do {                                                                         \
    __m128i keygend = _mm_setzero_si128();                                     \
    uint8_t s, t;                                                              \
    s = _mm_extract_epi8(m, 13);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    t = t ^ ROTL8(t, 1) ^ ROTL8(t, 2) ^ ROTL8(t, 3) ^ ROTL8(t, 4);             \
    keygend = _mm_insert_epi8(keygend, t ^ affine ^ rcon, 12);                 \
    s = _mm_extract_epi8(m, 14);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    t = t ^ ROTL8(t, 1) ^ ROTL8(t, 2) ^ ROTL8(t, 3) ^ ROTL8(t, 4);             \
    keygend = _mm_insert_epi8(keygend, t ^ affine, 13);                        \
    s = _mm_extract_epi8(m, 15);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    t = t ^ ROTL8(t, 1) ^ ROTL8(t, 2) ^ ROTL8(t, 3) ^ ROTL8(t, 4);             \
    keygend = _mm_insert_epi8(keygend, t ^ affine, 14);                        \
    s = _mm_extract_epi8(m, 12);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    t = t ^ ROTL8(t, 1) ^ ROTL8(t, 2) ^ ROTL8(t, 3) ^ ROTL8(t, 4);             \
    keygend = _mm_insert_epi8(keygend, t ^ affine, 15);                        \
    keygend = _mm_shuffle_epi32(keygend, _MM_SHUFFLE(3, 3, 3, 3));             \
    key = _mm_xor_si128(m, _mm_slli_si128(m, 4));                              \
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));                          \
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));                          \
    key = _mm_xor_si128(key, keygend);                                         \
  } while (0)

#define check_for_zeros_KS(m)                                                  \
  do {                                                                         \
    __m128i res = _mm_cmpeq_epi8(m, _mm_setzero_si128());                      \
    if (!_mm_test_all_zeros(res, _mm_set_epi32(-1, 0, 0, 0))) {                \
      return false;                                                            \
    }                                                                          \
  } while (0)

#define check_for_zeros(m)                                                     \
  do {                                                                         \
    __m128i res = _mm_cmpeq_epi8(m, _mm_setzero_si128());                      \
    if (!_mm_test_all_zeros(res, _mm_set1_epi64x(-1))) {                       \
      return false;                                                            \
    }                                                                          \
  } while (0)

// order is 13,14,15,12 due to rotword before sbox
#define check_for_zeros_and_save_KS(m, saved_sbox_state)                       \
  do {                                                                         \
    uint8_t s, t;                                                              \
    s = _mm_extract_epi8(m, 13);                                               \
    if (s == 0)                                                                \
      return false;                                                            \
    t = inverse(s);                                                            \
    saved_sbox_state.first.push_back(s);                                       \
    saved_sbox_state.second.push_back(t);                                      \
    s = _mm_extract_epi8(m, 14);                                               \
    if (s == 0)                                                                \
      return false;                                                            \
    t = inverse(s);                                                            \
    saved_sbox_state.first.push_back(s);                                       \
    saved_sbox_state.second.push_back(t);                                      \
    s = _mm_extract_epi8(m, 15);                                               \
    if (s == 0)                                                                \
      return false;                                                            \
    t = inverse(s);                                                            \
    saved_sbox_state.first.push_back(s);                                       \
    saved_sbox_state.second.push_back(t);                                      \
    s = _mm_extract_epi8(m, 12);                                               \
    if (s == 0)                                                                \
      return false;                                                            \
    t = inverse(s);                                                            \
    saved_sbox_state.first.push_back(s);                                       \
    saved_sbox_state.second.push_back(t);                                      \
  } while (0)

#define check_for_zeros_and_save(m, saved_sbox_state)                          \
  do {                                                                         \
    __m128i res = _mm_cmpeq_epi8(m, _mm_setzero_si128());                      \
    if (!_mm_test_all_zeros(res, _mm_set1_epi64x(-1))) {                       \
      return false;                                                            \
    }                                                                          \
    uint8_t buf[16];                                                           \
    _mm_storeu_si128((__m128i *)buf, m);                                       \
    for (int i = 0; i < 4; i++) {                                              \
      for (int j = 0; j < 4; j++) {                                            \
        saved_sbox_state.first.push_back(buf[j * 4 + i]);                      \
        saved_sbox_state.second.push_back(inverse(buf[j * 4 + i]));            \
      }                                                                        \
    }                                                                          \
  } while (0)

// order is 13,14,15,12 due to rotword before sbox
#define restore_t_shares_KS(m, s_shares, t_shares, party, sbox_index)          \
  do {                                                                         \
    uint8_t s, t;                                                              \
    s = _mm_extract_epi8(m, 13);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    s = inverse(t);                                                            \
    m = _mm_insert_epi8(m, s, 13);                                             \
    s = _mm_extract_epi8(m, 14);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    s = inverse(t);                                                            \
    m = _mm_insert_epi8(m, s, 14);                                             \
    s = _mm_extract_epi8(m, 15);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    s = inverse(t);                                                            \
    m = _mm_insert_epi8(m, s, 15);                                             \
    s = _mm_extract_epi8(m, 12);                                               \
    s_shares[party][sbox_index] = s;                                           \
    t = t_shares[party][sbox_index++];                                         \
    s = inverse(t);                                                            \
    m = _mm_insert_epi8(m, s, 12);                                             \
  } while (0)

#define restore_t_shares(m, s_shares, t_shares, party, sbox_index)             \
  do {                                                                         \
    uint8_t buf[16], buf2[16];                                                 \
    _mm_storeu_si128((__m128i *)buf, m);                                       \
    for (int i = 0; i < 4; i++) {                                              \
      for (int j = 0; j < 4; j++) {                                            \
        s_shares[party][sbox_index] = buf[j * 4 + i];                          \
        uint8_t t = t_shares[party][sbox_index++];                             \
        t = t ^ ROTL8(t, 1) ^ ROTL8(t, 2) ^ ROTL8(t, 3) ^ ROTL8(t, 4);         \
        if (party == 0)                                                        \
          t ^= AES_SBOX_AFFINE_CONST;                                          \
        buf2[((4 + j - i) % 4) * 4 + i] = t;                                   \
      }                                                                        \
    }                                                                          \
    m = _mm_loadu_si128((const __m128i *)buf2);                                \
  } while (0)

} // namespace

namespace AES128 {
static bool aes_128_aesni(const uint8_t *key, const uint8_t *plaintext,
                          uint8_t *ciphertext) {

  __m128i key_schedule[11];
  key_schedule[0] = _mm_loadu_si128((const __m128i *)key);
  check_for_zeros_KS(key_schedule[0]);
  key_schedule[1] = AES_128_key_exp(key_schedule[0], 0x01);
  check_for_zeros_KS(key_schedule[1]);
  key_schedule[2] = AES_128_key_exp(key_schedule[1], 0x02);
  check_for_zeros_KS(key_schedule[2]);
  key_schedule[3] = AES_128_key_exp(key_schedule[2], 0x04);
  check_for_zeros_KS(key_schedule[3]);
  key_schedule[4] = AES_128_key_exp(key_schedule[3], 0x08);
  check_for_zeros_KS(key_schedule[4]);
  key_schedule[5] = AES_128_key_exp(key_schedule[4], 0x10);
  check_for_zeros_KS(key_schedule[5]);
  key_schedule[6] = AES_128_key_exp(key_schedule[5], 0x20);
  check_for_zeros_KS(key_schedule[6]);
  key_schedule[7] = AES_128_key_exp(key_schedule[6], 0x40);
  check_for_zeros_KS(key_schedule[7]);
  key_schedule[8] = AES_128_key_exp(key_schedule[7], 0x80);
  check_for_zeros_KS(key_schedule[8]);
  key_schedule[9] = AES_128_key_exp(key_schedule[8], 0x1B);
  check_for_zeros_KS(key_schedule[9]);
  key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);

  __m128i m = _mm_loadu_si128((__m128i *)plaintext);
  m = _mm_xor_si128(m, key_schedule[0]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[1]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[2]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[3]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[4]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[5]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[6]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[7]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[8]);
  check_for_zeros(m);
  m = _mm_aesenc_si128(m, key_schedule[9]);
  check_for_zeros(m);
  m = _mm_aesenclast_si128(m, key_schedule[10]);

  _mm_storeu_si128((__m128i *)ciphertext, m);

  return true;
}

static bool aes_128_save_sbox_state_aesni(
    const uint8_t *key, const uint8_t *plaintext, uint8_t *ciphertext,
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> &saved_sbox_state) {
  __m128i key_schedule[11];
  key_schedule[0] = _mm_loadu_si128((const __m128i *)key);
  check_for_zeros_and_save_KS(key_schedule[0], saved_sbox_state);
  key_schedule[1] = AES_128_key_exp(key_schedule[0], 0x01);
  check_for_zeros_and_save_KS(key_schedule[1], saved_sbox_state);
  key_schedule[2] = AES_128_key_exp(key_schedule[1], 0x02);
  check_for_zeros_and_save_KS(key_schedule[2], saved_sbox_state);
  key_schedule[3] = AES_128_key_exp(key_schedule[2], 0x04);
  check_for_zeros_and_save_KS(key_schedule[3], saved_sbox_state);
  key_schedule[4] = AES_128_key_exp(key_schedule[3], 0x08);
  check_for_zeros_and_save_KS(key_schedule[4], saved_sbox_state);
  key_schedule[5] = AES_128_key_exp(key_schedule[4], 0x10);
  check_for_zeros_and_save_KS(key_schedule[5], saved_sbox_state);
  key_schedule[6] = AES_128_key_exp(key_schedule[5], 0x20);
  check_for_zeros_and_save_KS(key_schedule[6], saved_sbox_state);
  key_schedule[7] = AES_128_key_exp(key_schedule[6], 0x40);
  check_for_zeros_and_save_KS(key_schedule[7], saved_sbox_state);
  key_schedule[8] = AES_128_key_exp(key_schedule[7], 0x80);
  check_for_zeros_and_save_KS(key_schedule[8], saved_sbox_state);
  key_schedule[9] = AES_128_key_exp(key_schedule[8], 0x1B);
  check_for_zeros_and_save_KS(key_schedule[9], saved_sbox_state);
  key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);

  __m128i m = _mm_loadu_si128((__m128i *)plaintext);
  m = _mm_xor_si128(m, key_schedule[0]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[1]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[2]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[3]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[4]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[5]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[6]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[7]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[8]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenc_si128(m, key_schedule[9]);
  check_for_zeros_and_save(m, saved_sbox_state);
  m = _mm_aesenclast_si128(m, key_schedule[10]);

  _mm_storeu_si128((__m128i *)ciphertext, m);

  return true;
}

bool aes_128(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintext_in,
             std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES128::KEY_SIZE);
  assert(plaintext_in.size() == AES128::BLOCK_SIZE * AES128::NUM_BLOCKS);
  ciphertext_out.resize(AES128::BLOCK_SIZE * AES128::NUM_BLOCKS);
  return aes_128_aesni(key_in.data(), plaintext_in.data(),
                       ciphertext_out.data());
}

std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
aes_128_with_sbox_output(const std::vector<uint8_t> &key_in,
                         const std::vector<uint8_t> &plaintext_in,
                         std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES128::KEY_SIZE);
  assert(plaintext_in.size() == AES128::BLOCK_SIZE * AES128::NUM_BLOCKS);
  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> result;
  result.first.reserve(AES128::NUM_SBOXES);
  result.second.reserve(AES128::NUM_SBOXES);
  ciphertext_out.resize(AES128::BLOCK_SIZE * AES128::NUM_BLOCKS);
  bool ret = aes_128_save_sbox_state_aesni(key_in.data(), plaintext_in.data(),
                                           ciphertext_out.data(), result);
  (void)ret;
  assert(ret);
  return result;
}

void aes_128_s_shares_old(const std::vector<gsl::span<uint8_t>> &key_in,
                          const std::vector<gsl::span<uint8_t>> &t_shares,
                          const std::vector<uint8_t> &plaintext_in,
                          std::vector<gsl::span<uint8_t>> &ciphertext_out,
                          std::vector<gsl::span<uint8_t>> &s_shares) {

  typedef std::array<std::array<uint8_t, 44>, 4> expanded_key_t;
  typedef std::array<std::array<uint8_t, 4>, 4> state_t;
  typedef std::array<uint8_t, 4> temp_t;
  int num_parties = key_in.size();
  std::vector<expanded_key_t> expanded(num_parties);
  std::vector<state_t> state(num_parties);
  std::vector<state_t> newstate(num_parties);
  uint8_t roundconstant;
  int i;
  int j;
  int r;
  int party;
  int sbox_index = 0;
  int first_party = 0;

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      for (party = 0; party < num_parties; party++)
        expanded[party][i][j] = key_in[party][j * 4 + i];

  roundconstant = 1;
  for (j = 4; j < 44; ++j) {
    std::vector<temp_t> temp(num_parties);
    if (j % 4)
      for (i = 0; i < 4; ++i)
        for (party = 0; party < num_parties; party++)
          temp[party][i] = expanded[party][i][j - 1];
    else {
      for (i = 0; i < 4; ++i) {
        for (party = 0; party < num_parties; party++) {
          unsigned char s = expanded[party][(i + 1) % 4][j - 1];
          s_shares[party][sbox_index] = s;
          temp[party][i] = bytesub_restore(t_shares[party][sbox_index]);
        }
        temp[0][i] ^= AES_SBOX_AFFINE_CONST;
        sbox_index++;
      }
      temp[first_party][0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      for (party = 0; party < num_parties; party++) {
        expanded[party][i][j] = temp[party][i] ^ expanded[party][i][j - 4];
      }
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i) {
      state[first_party][i][j] =
          plaintext_in[j * 4 + i] ^ expanded[first_party][i][j];
      for (party = 1; party < num_parties; party++) {
        state[party][i][j] = expanded[party][i][j];
      }
    }

  for (r = 0; r < 10; ++r) {
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
        for (party = 0; party < num_parties; party++) {
          s_shares[party][sbox_index] = state[party][i][j];
          newstate[party][i][j] = bytesub_restore(t_shares[party][sbox_index]);
        }
        newstate[0][i][j] ^= AES_SBOX_AFFINE_CONST;
        sbox_index++;
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        for (party = 0; party < num_parties; party++) {
          state[party][i][j] = newstate[party][i][(j + i) % 4];
        }
    if (r < 9)
      for (j = 0; j < 4; ++j) {
        for (party = 0; party < num_parties; party++) {
          unsigned char a0 = state[party][0][j];
          unsigned char a1 = state[party][1][j];
          unsigned char a2 = state[party][2][j];
          unsigned char a3 = state[party][3][j];
          state[party][0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
          state[party][1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
          state[party][2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
          state[party][3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
        }
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        for (party = 0; party < num_parties; party++) {
          state[party][i][j] ^= expanded[party][i][r * 4 + 4 + j];
        }
  }

  for (party = 0; party < num_parties; party++) {
    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i)
        ciphertext_out[party][j * 4 + i] = state[party][i][j];
  }
}
void aes_128_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
                      const std::vector<gsl::span<uint8_t>> &t_shares,
                      const std::vector<uint8_t> &plaintext_in,
                      std::vector<gsl::span<uint8_t>> &ciphertext_out,
                      std::vector<gsl::span<uint8_t>> &s_shares) {

  typedef std::array<__m128i, 11> expanded_key_t;
  int num_parties = key_in.size();
  std::vector<expanded_key_t> key_schedule(num_parties);
  std::vector<__m128i> state(num_parties);
  int party = 0;
  int sbox_index = 0;
  // first party do normal sbox + rcon
  {
    key_schedule[party][0] =
        _mm_loadu_si128((const __m128i *)key_in[party].data());
    AES_128_key_exp_restore(key_schedule[party][1], key_schedule[party][0],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x1);
    AES_128_key_exp_restore(key_schedule[party][2], key_schedule[party][1],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x2);
    AES_128_key_exp_restore(key_schedule[party][3], key_schedule[party][2],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x4);
    AES_128_key_exp_restore(key_schedule[party][4], key_schedule[party][3],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x8);
    AES_128_key_exp_restore(key_schedule[party][5], key_schedule[party][4],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x10);
    AES_128_key_exp_restore(key_schedule[party][6], key_schedule[party][5],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x20);
    AES_128_key_exp_restore(key_schedule[party][7], key_schedule[party][6],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x40);
    AES_128_key_exp_restore(key_schedule[party][8], key_schedule[party][7],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x80);
    AES_128_key_exp_restore(key_schedule[party][9], key_schedule[party][8],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x1B);
    AES_128_key_exp_restore(key_schedule[party][10], key_schedule[party][9],
                            s_shares, t_shares, party, sbox_index,
                            AES_SBOX_AFFINE_CONST, 0x36);
  }
  // other parties do not add rcon or affine part of sbox
  for (party = 1; party < num_parties; party++) {
    sbox_index = 0;
    key_schedule[party][0] =
        _mm_loadu_si128((const __m128i *)key_in[party].data());
    AES_128_key_exp_restore(key_schedule[party][1], key_schedule[party][0],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][2], key_schedule[party][1],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][3], key_schedule[party][2],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][4], key_schedule[party][3],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][5], key_schedule[party][4],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][6], key_schedule[party][5],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][7], key_schedule[party][6],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][8], key_schedule[party][7],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][9], key_schedule[party][8],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
    AES_128_key_exp_restore(key_schedule[party][10], key_schedule[party][9],
                            s_shares, t_shares, party, sbox_index, 0x0, 0x0);
  }
  // fix up the additionall affine part in all parties but the first one
  // that arises from actually using aes instructions
  // the affine part of the sbox should only happen at party 0, so
  // we cancel out the additional part for the other parties after MixCol
  for (party = 0; party < num_parties; party++) {
    sbox_index = 40;
    if (party == 0)
      state[party] = _mm_loadu_si128((__m128i *)plaintext_in.data());
    else
      state[party] = _mm_setzero_si128();
    state[party] = _mm_xor_si128(state[party], key_schedule[party][0]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][1]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][2]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][3]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][4]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][5]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][6]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][7]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][8]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] =
        _mm_aesimc_si128(_mm_aesimc_si128(_mm_aesimc_si128(state[party])));
    state[party] = _mm_xor_si128(state[party], key_schedule[party][9]);
    restore_t_shares(state[party], s_shares, t_shares, party, sbox_index);
    state[party] = _mm_xor_si128(state[party], key_schedule[party][10]);

    _mm_storeu_si128((__m128i *)ciphertext_out[party].data(), state[party]);
  }
}

} // namespace AES128

namespace AES192 {

static bool aes_192_old(const uint8_t *key, const uint8_t *plaintext,
                        uint8_t *ciphertext) {
  unsigned char expanded[4][52];
  unsigned char state[4][4];
  unsigned char newstate[4][4];
  unsigned char roundconstant;
  int i;
  int j;
  int r;

  for (j = 0; j < 6; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key[j * 4 + i];

  roundconstant = 1;
  for (j = 6; j < 52; ++j) {
    unsigned char temp[4];
    if (j % 6)
      for (i = 0; i < 4; ++i)
        temp[i] = expanded[i][j - 1];
    else {
      for (i = 0; i < 4; ++i) {
        unsigned char s = expanded[(i + 1) % 4][j - 1];
        if (s == 0)
          return false;
        temp[i] = bytesub(s);
      }
      temp[0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      expanded[i][j] = temp[i] ^ expanded[i][j - 6];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      state[i][j] = plaintext[j * 4 + i] ^ expanded[i][j];

  for (r = 0; r < 12; ++r) {
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
        if (state[i][j] == 0)
          return false;
        newstate[i][j] = bytesub(state[i][j]);
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] = newstate[i][(j + i) % 4];
    if (r < 11)
      for (j = 0; j < 4; ++j) {
        unsigned char a0 = state[0][j];
        unsigned char a1 = state[1][j];
        unsigned char a2 = state[2][j];
        unsigned char a3 = state[3][j];
        state[0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
        state[1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
        state[2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
        state[3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] ^= expanded[i][r * 4 + 4 + j];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      ciphertext[j * 4 + i] = state[i][j];

  return true;
}

static bool aes_192_save_sbox_state(
    const uint8_t *key, const uint8_t *plaintext, uint8_t *ciphertext,
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> &saved_sbox_state) {
  unsigned char expanded[4][52];
  unsigned char state[4][4];
  unsigned char newstate[4][4];
  unsigned char roundconstant;
  int i;
  int j;
  int r;

  for (j = 0; j < 6; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key[j * 4 + i];

  roundconstant = 1;
  for (j = 6; j < 52; ++j) {
    unsigned char temp[4];
    if (j % 6)
      for (i = 0; i < 4; ++i)
        temp[i] = expanded[i][j - 1];
    else {
      for (i = 0; i < 4; ++i) {
        unsigned char s = expanded[(i + 1) % 4][j - 1];
        if (s == 0)
          return false;
        std::pair<uint8_t, uint8_t> sbox_state;
        temp[i] = bytesub_save(s, sbox_state);
        saved_sbox_state.first.push_back(sbox_state.first);
        saved_sbox_state.second.push_back(sbox_state.second);
      }
      temp[0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      expanded[i][j] = temp[i] ^ expanded[i][j - 6];
  }

  for (size_t k = 0; k < AES192::NUM_BLOCKS; k++) {
    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i)
        state[i][j] =
            plaintext[k * AES192::BLOCK_SIZE + j * 4 + i] ^ expanded[i][j];

    for (r = 0; r < 12; ++r) {
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j) {
          if (state[i][j] == 0)
            return false;
          std::pair<uint8_t, uint8_t> sbox_state;
          newstate[i][j] = bytesub_save(state[i][j], sbox_state);
          saved_sbox_state.first.push_back(sbox_state.first);
          saved_sbox_state.second.push_back(sbox_state.second);
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          state[i][j] = newstate[i][(j + i) % 4];
      if (r < 11)
        for (j = 0; j < 4; ++j) {
          unsigned char a0 = state[0][j];
          unsigned char a1 = state[1][j];
          unsigned char a2 = state[2][j];
          unsigned char a3 = state[3][j];
          state[0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
          state[1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
          state[2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
          state[3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          state[i][j] ^= expanded[i][r * 4 + 4 + j];
    }

    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i)
        ciphertext[k * AES192::BLOCK_SIZE + j * 4 + i] = state[i][j];
  }

  return true;
}

bool aes_192(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintext_in,
             std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES192::KEY_SIZE);
  assert(plaintext_in.size() == AES192::BLOCK_SIZE * AES192::NUM_BLOCKS);
  ciphertext_out.resize(AES192::BLOCK_SIZE * AES192::NUM_BLOCKS);
  bool ok =
      aes_192_old(key_in.data(), plaintext_in.data(), ciphertext_out.data());
  ok = aes_192_old(key_in.data(), plaintext_in.data() + AES192::BLOCK_SIZE,
                   ciphertext_out.data() + AES192::BLOCK_SIZE) &&
       ok;
  return ok;
}
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
aes_192_with_sbox_output(const std::vector<uint8_t> &key_in,
                         const std::vector<uint8_t> &plaintext_in,
                         std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES192::KEY_SIZE);
  assert(plaintext_in.size() == AES192::BLOCK_SIZE * AES192::NUM_BLOCKS);
  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> result;
  result.first.reserve(AES192::NUM_SBOXES);
  result.second.reserve(AES192::NUM_SBOXES);
  ciphertext_out.resize(AES192::BLOCK_SIZE * AES192::NUM_BLOCKS);
  bool ret = aes_192_save_sbox_state(key_in.data(), plaintext_in.data(),
                                     ciphertext_out.data(), result);
  (void)ret;
  assert(ret);
  return result;
}

void aes_192_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
                      const std::vector<gsl::span<uint8_t>> &t_shares,
                      const std::vector<uint8_t> &plaintext_in,
                      std::vector<gsl::span<uint8_t>> &ciphertext_out,
                      std::vector<gsl::span<uint8_t>> &s_shares) {

  typedef std::array<std::array<uint8_t, 52>, 4> expanded_key_t;
  typedef std::array<std::array<uint8_t, 4>, 4> state_t;
  typedef std::array<uint8_t, 4> temp_t;
  int num_parties = key_in.size();
  std::vector<expanded_key_t> expanded(num_parties);
  std::vector<state_t> state(num_parties);
  std::vector<state_t> newstate(num_parties);
  uint8_t roundconstant;
  int i;
  int j;
  int r;
  int party;
  int sbox_index = 0;
  int first_party = 0;

  for (j = 0; j < 6; ++j)
    for (i = 0; i < 4; ++i)
      for (party = 0; party < num_parties; party++)
        expanded[party][i][j] = key_in[party][j * 4 + i];

  roundconstant = 1;
  for (j = 6; j < 52; ++j) {
    std::vector<temp_t> temp(num_parties);
    if (j % 6)
      for (i = 0; i < 4; ++i)
        for (party = 0; party < num_parties; party++)
          temp[party][i] = expanded[party][i][j - 1];
    else {
      for (i = 0; i < 4; ++i) {
        for (party = 0; party < num_parties; party++) {
          unsigned char s = expanded[party][(i + 1) % 4][j - 1];
          s_shares[party][sbox_index] = s;
          temp[party][i] = bytesub_restore(t_shares[party][sbox_index]);
        }
        temp[0][i] ^= AES_SBOX_AFFINE_CONST;
        sbox_index++;
      }
      temp[first_party][0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      for (party = 0; party < num_parties; party++) {
        expanded[party][i][j] = temp[party][i] ^ expanded[party][i][j - 6];
      }
  }

  for (size_t k = 0; k < AES192::NUM_BLOCKS; k++) {
    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i) {
        state[first_party][i][j] =
            plaintext_in[k * AES192::BLOCK_SIZE + j * 4 + i] ^
            expanded[first_party][i][j];
        for (party = 1; party < num_parties; party++) {
          state[party][i][j] = expanded[party][i][j];
        }
      }

    for (r = 0; r < 12; ++r) {
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j) {
          for (party = 0; party < num_parties; party++) {
            s_shares[party][sbox_index] = state[party][i][j];
            newstate[party][i][j] =
                bytesub_restore(t_shares[party][sbox_index]);
          }
          newstate[0][i][j] ^= AES_SBOX_AFFINE_CONST;
          sbox_index++;
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          for (party = 0; party < num_parties; party++) {
            state[party][i][j] = newstate[party][i][(j + i) % 4];
          }
      if (r < 11)
        for (j = 0; j < 4; ++j) {
          for (party = 0; party < num_parties; party++) {
            unsigned char a0 = state[party][0][j];
            unsigned char a1 = state[party][1][j];
            unsigned char a2 = state[party][2][j];
            unsigned char a3 = state[party][3][j];
            state[party][0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
            state[party][1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
            state[party][2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
            state[party][3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
          }
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          for (party = 0; party < num_parties; party++) {
            state[party][i][j] ^= expanded[party][i][r * 4 + 4 + j];
          }
    }

    for (party = 0; party < num_parties; party++) {
      for (j = 0; j < 4; ++j)
        for (i = 0; i < 4; ++i)
          ciphertext_out[party][k * AES192::BLOCK_SIZE + j * 4 + i] =
              state[party][i][j];
    }
  }
}
} // namespace AES192

namespace AES256 {
bool aes_256_old(const unsigned char *key, const unsigned char *plaintext,
                 unsigned char *ciphertext) {
  unsigned char expanded[4][60];
  unsigned char state[4][4];
  unsigned char newstate[4][4];
  unsigned char roundconstant;
  int i;
  int j;
  int r;

  for (j = 0; j < 8; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key[j * 4 + i];

  roundconstant = 1;
  for (j = 8; j < 60; ++j) {
    unsigned char temp[4];
    if (j % 4)
      for (i = 0; i < 4; ++i)
        temp[i] = expanded[i][j - 1];
    else if (j % 8)
      for (i = 0; i < 4; ++i) {
        if (expanded[i][j - 1] == 0)
          return false;
        temp[i] = bytesub(expanded[i][j - 1]);
      }
    else {
      for (i = 0; i < 4; ++i) {
        if (expanded[(i + 1) % 4][j - 1] == 0)
          return false;
        temp[i] = bytesub(expanded[(i + 1) % 4][j - 1]);
      }
      temp[0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      expanded[i][j] = temp[i] ^ expanded[i][j - 8];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      state[i][j] = plaintext[j * 4 + i] ^ expanded[i][j];

  for (r = 0; r < 14; ++r) {
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
        if (state[i][j] == 0)
          return false;
        newstate[i][j] = bytesub(state[i][j]);
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] = newstate[i][(j + i) % 4];
    if (r < 13)
      for (j = 0; j < 4; ++j) {
        unsigned char a0 = state[0][j];
        unsigned char a1 = state[1][j];
        unsigned char a2 = state[2][j];
        unsigned char a3 = state[3][j];
        state[0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
        state[1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
        state[2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
        state[3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] ^= expanded[i][r * 4 + 4 + j];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      ciphertext[j * 4 + i] = state[i][j];

  return true;
}

bool aes_256_save_sbox_state(
    const unsigned char *key, const unsigned char *plaintext,
    unsigned char *ciphertext,
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> &saved_sbox_state) {
  unsigned char expanded[4][60];
  unsigned char state[4][4];
  unsigned char newstate[4][4];
  unsigned char roundconstant;
  int i;
  int j;
  int r;

  for (j = 0; j < 8; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key[j * 4 + i];

  roundconstant = 1;
  for (j = 8; j < 60; ++j) {
    unsigned char temp[4];
    if (j % 4)
      for (i = 0; i < 4; ++i)
        temp[i] = expanded[i][j - 1];
    else if (j % 8)
      for (i = 0; i < 4; ++i) {
        if (expanded[i][j - 1] == 0)
          return false;
        std::pair<uint8_t, uint8_t> sbox_state;
        temp[i] = bytesub_save(expanded[i][j - 1], sbox_state);
        saved_sbox_state.first.push_back(sbox_state.first);
        saved_sbox_state.second.push_back(sbox_state.second);
      }
    else {
      for (i = 0; i < 4; ++i) {
        if (expanded[(i + 1) % 4][j - 1] == 0)
          return false;
        std::pair<uint8_t, uint8_t> sbox_state;
        temp[i] = bytesub_save(expanded[(i + 1) % 4][j - 1], sbox_state);
        saved_sbox_state.first.push_back(sbox_state.first);
        saved_sbox_state.second.push_back(sbox_state.second);
      }
      temp[0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      expanded[i][j] = temp[i] ^ expanded[i][j - 8];
  }

  for (size_t k = 0; k < AES256::NUM_BLOCKS; k++) {
    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i)
        state[i][j] =
            plaintext[k * AES256::BLOCK_SIZE + j * 4 + i] ^ expanded[i][j];

    for (r = 0; r < 14; ++r) {
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j) {
          if (state[i][j] == 0)
            return false;
          std::pair<uint8_t, uint8_t> sbox_state;
          newstate[i][j] = bytesub_save(state[i][j], sbox_state);
          saved_sbox_state.first.push_back(sbox_state.first);
          saved_sbox_state.second.push_back(sbox_state.second);
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          state[i][j] = newstate[i][(j + i) % 4];
      if (r < 13)
        for (j = 0; j < 4; ++j) {
          unsigned char a0 = state[0][j];
          unsigned char a1 = state[1][j];
          unsigned char a2 = state[2][j];
          unsigned char a3 = state[3][j];
          state[0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
          state[1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
          state[2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
          state[3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          state[i][j] ^= expanded[i][r * 4 + 4 + j];
    }

    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i)
        ciphertext[k * AES256::BLOCK_SIZE + j * 4 + i] = state[i][j];
  }

  return true;
}
bool aes_256(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintext_in,
             std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES256::KEY_SIZE);
  assert(plaintext_in.size() == AES256::BLOCK_SIZE * AES256::NUM_BLOCKS);
  ciphertext_out.resize(AES128::BLOCK_SIZE * AES256::NUM_BLOCKS);
  bool ok =
      aes_256_old(key_in.data(), plaintext_in.data(), ciphertext_out.data());
  ok = aes_256_old(key_in.data(), plaintext_in.data() + AES256::BLOCK_SIZE,
                   ciphertext_out.data() + AES256::BLOCK_SIZE) &&
       ok;
  return ok;
}
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
aes_256_with_sbox_output(const std::vector<uint8_t> &key_in,
                         const std::vector<uint8_t> &plaintext_in,
                         std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES256::KEY_SIZE);
  assert(plaintext_in.size() == AES256::BLOCK_SIZE * AES256::NUM_BLOCKS);
  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> result;
  result.first.reserve(AES256::NUM_SBOXES);
  result.second.reserve(AES256::NUM_SBOXES);
  ciphertext_out.resize(AES256::BLOCK_SIZE * AES256::NUM_BLOCKS);
  bool ret = aes_256_save_sbox_state(key_in.data(), plaintext_in.data(),
                                     ciphertext_out.data(), result);

  // remove duplicated keyschedule from saved sbox values
  // keyschedule = 52 sboxes, rounds = 224 delete from
  (void)ret;
  assert(ret);
  return result;
}

void aes_256_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
                      const std::vector<gsl::span<uint8_t>> &t_shares,
                      const std::vector<uint8_t> &plaintext_in,
                      std::vector<gsl::span<uint8_t>> &ciphertext_out,
                      std::vector<gsl::span<uint8_t>> &s_shares) {

  typedef std::array<std::array<uint8_t, 60>, 4> expanded_key_t;
  typedef std::array<std::array<uint8_t, 4>, 4> state_t;
  typedef std::array<uint8_t, 4> temp_t;
  int num_parties = key_in.size();
  std::vector<expanded_key_t> expanded(num_parties);
  std::vector<state_t> state(num_parties);
  std::vector<state_t> newstate(num_parties);
  uint8_t roundconstant;
  int i;
  int j;
  int r;
  int party;
  int sbox_index = 0;
  int first_party = 0;

  for (j = 0; j < 8; ++j)
    for (i = 0; i < 4; ++i)
      for (party = 0; party < num_parties; party++)
        expanded[party][i][j] = key_in[party][j * 4 + i];

  roundconstant = 1;
  for (j = 8; j < 60; ++j) {
    std::vector<temp_t> temp(num_parties);
    if (j % 4)
      for (i = 0; i < 4; ++i)
        for (party = 0; party < num_parties; party++)
          temp[party][i] = expanded[party][i][j - 1];
    else if (j % 8)
      for (i = 0; i < 4; ++i) {
        for (party = 0; party < num_parties; party++) {
          unsigned char s = expanded[party][i][j - 1];
          s_shares[party][sbox_index] = s;
          temp[party][i] = bytesub_restore(t_shares[party][sbox_index]);
        }
        temp[0][i] ^= AES_SBOX_AFFINE_CONST;
        sbox_index++;
      }
    else {
      for (i = 0; i < 4; ++i) {
        for (party = 0; party < num_parties; party++) {
          unsigned char s = expanded[party][(i + 1) % 4][j - 1];
          s_shares[party][sbox_index] = s;
          temp[party][i] = bytesub_restore(t_shares[party][sbox_index]);
        }
        temp[0][i] ^= AES_SBOX_AFFINE_CONST;
        sbox_index++;
      }
      temp[first_party][0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      for (party = 0; party < num_parties; party++) {
        expanded[party][i][j] = temp[party][i] ^ expanded[party][i][j - 8];
      }
  }

  for (size_t k = 0; k < AES256::NUM_BLOCKS; k++) {
    for (j = 0; j < 4; ++j)
      for (i = 0; i < 4; ++i) {
        state[first_party][i][j] =
            plaintext_in[k * AES256::BLOCK_SIZE + j * 4 + i] ^
            expanded[first_party][i][j];
        for (party = 1; party < num_parties; party++) {
          state[party][i][j] = expanded[party][i][j];
        }
      }

    for (r = 0; r < 14; ++r) {
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j) {
          for (party = 0; party < num_parties; party++) {
            s_shares[party][sbox_index] = state[party][i][j];
            newstate[party][i][j] =
                bytesub_restore(t_shares[party][sbox_index]);
          }
          newstate[0][i][j] ^= AES_SBOX_AFFINE_CONST;
          sbox_index++;
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          for (party = 0; party < num_parties; party++) {
            state[party][i][j] = newstate[party][i][(j + i) % 4];
          }
      if (r < 13)
        for (j = 0; j < 4; ++j) {
          for (party = 0; party < num_parties; party++) {
            unsigned char a0 = state[party][0][j];
            unsigned char a1 = state[party][1][j];
            unsigned char a2 = state[party][2][j];
            unsigned char a3 = state[party][3][j];
            state[party][0][j] = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
            state[party][1][j] = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
            state[party][2][j] = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
            state[party][3][j] = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
          }
        }
      for (i = 0; i < 4; ++i)
        for (j = 0; j < 4; ++j)
          for (party = 0; party < num_parties; party++) {
            state[party][i][j] ^= expanded[party][i][r * 4 + 4 + j];
          }
    }
    for (party = 0; party < num_parties; party++) {
      for (j = 0; j < 4; ++j)
        for (i = 0; i < 4; ++i)
          ciphertext_out[party][k * AES256::BLOCK_SIZE + j * 4 + i] =
              state[party][i][j];
    }
  }
}
} // namespace AES256