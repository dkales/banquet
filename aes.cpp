#include "aes.h"
#include <cassert>

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
} // namespace

namespace AES128 {
static bool aes_128_old(const uint8_t *key, const uint8_t *plaintext,
                        uint8_t *ciphertext) {
  unsigned char expanded[4][44];
  unsigned char state[4][4];
  unsigned char newstate[4][4];
  unsigned char roundconstant;
  int i;
  int j;
  int r;

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key[j * 4 + i];

  roundconstant = 1;
  for (j = 4; j < 44; ++j) {
    unsigned char temp[4];
    if (j % 4)
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
      expanded[i][j] = temp[i] ^ expanded[i][j - 4];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      state[i][j] = plaintext[j * 4 + i] ^ expanded[i][j];

  for (r = 0; r < 10; ++r) {
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
        if (state[i][j] == 0)
          return false;
        newstate[i][j] = bytesub(state[i][j]);
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] = newstate[i][(j + i) % 4];
    if (r < 9)
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

static bool aes_128_save_sbox_state(
    const uint8_t *key, const uint8_t *plaintext, uint8_t *ciphertext,
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> &saved_sbox_state) {
  unsigned char expanded[4][44];
  unsigned char state[4][4];
  unsigned char newstate[4][4];
  unsigned char roundconstant;
  int i;
  int j;
  int r;

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key[j * 4 + i];

  roundconstant = 1;
  for (j = 4; j < 44; ++j) {
    unsigned char temp[4];
    if (j % 4)
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
      expanded[i][j] = temp[i] ^ expanded[i][j - 4];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      state[i][j] = plaintext[j * 4 + i] ^ expanded[i][j];

  for (r = 0; r < 10; ++r) {
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
    if (r < 9)
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

bool aes_128(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintext_in,
             std::vector<uint8_t> &ciphertext_out) {
  assert(key_in.size() == AES128::KEY_SIZE);
  assert(plaintext_in.size() == AES128::BLOCK_SIZE * AES128::NUM_BLOCKS);
  ciphertext_out.resize(AES128::BLOCK_SIZE * AES128::NUM_BLOCKS);
  return aes_128_old(key_in.data(), plaintext_in.data(), ciphertext_out.data());
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
  bool ret = aes_128_save_sbox_state(key_in.data(), plaintext_in.data(),
                                     ciphertext_out.data(), result);
  (void)ret;
  assert(ret);
  return result;
}

void aes_128_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
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