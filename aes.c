#include "aes.h"
#include "io.h"
#include <string.h>

int aes_128(const mzd_local_t* key_in, const mzd_local_t* plaintext_in,
            mzd_local_t* ciphertext_out) {
  uint8_t plaintext[16];
  uint8_t ciphertext[16];
  uint8_t key[16];

  mzd_to_char_array(key, key_in, 16);
  mzd_to_char_array(plaintext, plaintext_in, 16);

  int ret = aes_128_old(key, plaintext, ciphertext);
  mzd_from_char_array(ciphertext_out, ciphertext, 16);
  return ret;
}

static unsigned char multiply(unsigned int a, unsigned int b) {
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

static unsigned char multiply2(unsigned int c, unsigned int d) {
  unsigned char f[8];
  unsigned char g[8];
  unsigned char h[15];
  unsigned char result;
  int i;
  int j;

  for (i = 0; i < 8; ++i)
    f[i] = 1 & (c >> i);
  for (i = 0; i < 8; ++i)
    g[i] = 1 & (d >> i);
  for (i = 0; i < 15; ++i)
    h[i] = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      h[i + j] ^= f[i] & g[j];

  for (i = 6; i >= 0; --i) {
    h[i + 0] ^= h[i + 8];
    h[i + 1] ^= h[i + 8];
    h[i + 3] ^= h[i + 8];
    h[i + 4] ^= h[i + 8];
    h[i + 8] ^= h[i + 8];
  }

  result = 0;
  for (i = 0; i < 8; ++i)
    result |= h[i] << i;
  return result;
}

static unsigned char square(unsigned char c) {
  return multiply(c, c);
}

static unsigned char xtime(unsigned char c) {
  return multiply(c, 2);
}

static unsigned char bytesub(unsigned char c) {
  unsigned char c3   = multiply(square(c), c);
  unsigned char c7   = multiply(square(c3), c);
  unsigned char c63  = multiply(square(square(square(c7))), c7);
  unsigned char c127 = multiply(square(c63), c);
  unsigned char c254 = square(c127);
  unsigned char f[8];
  unsigned char h[8];
  unsigned char result = c254;
  int i;
  
  for (i = 0; i < 8; ++i)
    f[i] = 1 & (c254 >> i);
  h[0]   = f[0] ^ f[4] ^ f[5] ^ f[6] ^ f[7] ^ 1;
  h[1]   = f[1] ^ f[5] ^ f[6] ^ f[7] ^ f[0] ^ 1;
  h[2]   = f[2] ^ f[6] ^ f[7] ^ f[0] ^ f[1];
  h[3]   = f[3] ^ f[7] ^ f[0] ^ f[1] ^ f[2];
  h[4]   = f[4] ^ f[0] ^ f[1] ^ f[2] ^ f[3];
  h[5]   = f[5] ^ f[1] ^ f[2] ^ f[3] ^ f[4] ^ 1;
  h[6]   = f[6] ^ f[2] ^ f[3] ^ f[4] ^ f[5] ^ 1;
  h[7]   = f[7] ^ f[3] ^ f[4] ^ f[5] ^ f[6];
  result = 0;
  for (i = 0; i < 8; ++i)
    result |= h[i] << i;
  //printf("%u->%u\n", c, result);
  return result;
}

int aes_128_old(const uint8_t* key, const uint8_t* plaintext, uint8_t* ciphertext) {
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
        if (expanded[(i + 1) % 4][j - 1] == 0)
          return 1;
        temp[i] = bytesub(expanded[(i + 1) % 4][j - 1]);
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
          return 1;
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
        state[0][j]      = xtime(a0 ^ a1) ^ a1 ^ a2 ^ a3;
        state[1][j]      = xtime(a1 ^ a2) ^ a2 ^ a3 ^ a0;
        state[2][j]      = xtime(a2 ^ a3) ^ a3 ^ a0 ^ a1;
        state[3][j]      = xtime(a3 ^ a0) ^ a0 ^ a1 ^ a2;
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] ^= expanded[i][r * 4 + 4 + j];
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      ciphertext[j * 4 + i] = state[i][j];

  return 0;
}

void aux_aes_128(const mzd_local_t* key, randomTape_t* tapes) {
  (void)key;
  const size_t NUM_AES_128_SBOXES = 200;
  const size_t BUFFER             = 10;

  // skip first 16 bytes for sharing of key
  tapes->pos = 16;
  for (size_t i = 0; i < NUM_AES_128_SBOXES + BUFFER; i++) {
    uint8_t a, b, c;
    uint32_t abc = 0;
    for (size_t j = 0; j < 64; j++) {
      // 4*i + 0 is the random byte r used in inversion, this is ok to be fully random, so we dont
      // fix it
      //a ^= tapes->tape[j][16 + 4 * i + 1];
      //b ^= tapes->tape[j][16 + 4 * i + 2];
      //c ^= tapes->tape[j][16 + 4 * i + 3];
      abc ^= ((uint32_t*)(tapes->tape[j]))[4 + i];
    }
    a             = (abc >> 8);
    b             = (abc >> 16);
    c             = (abc >> 24);
    uint8_t c_err = multiply(a, b) ^ c;
    tapes->tape[63][16 + 4 * i + 3] ^= c_err;
    tapes->aux_bits[i] = tapes->tape[63][16 + 4 * i + 3];
  }
  tapes->pos = 0;
}

static inline aes_byte_share_t getByteFromTapes(randomTape_t* tapes) {
  aes_byte_share_t ret;
  for (size_t i = 0; i < 64; i++) {
    ret.byte[i] = tapes->tape[i][tapes->pos];
  }
  tapes->pos++;
  return ret;
}

static inline aes_byte_share_t subShare(aes_byte_share_t a, aes_byte_share_t b) {
  aes_byte_share_t ret;
  for (size_t i = 0; i < 64; i++) {
    ret.byte[i] = a.byte[i] ^ b.byte[i];
  }
  return ret;
}

static inline uint8_t broadcastShare(aes_byte_share_t a, msgs_t* msgs) {
  uint8_t ret = 0;
  if (msgs->unopened >= 0) {
    a.byte[msgs->unopened] = msgs->msgs[msgs->unopened][msgs->pos / 8];
  }
  for (size_t i = 0; i < 64; i++) {
    msgs->msgs[i][msgs->pos/8] = a.byte[i];
    ret ^= a.byte[i];
  }
  msgs->pos+=8;
  return ret;
}

static inline aes_byte_share_t computeMul(aes_byte_share_t a, aes_byte_share_t b,
                                          aes_byte_share_t c, uint8_t d, uint8_t e) {
  aes_byte_share_t ret;
  uint8_t de = multiply(d, e);
  for (size_t i = 0; i < 64; i++) {
    ret.byte[i] = c.byte[i] ^ multiply(d, b.byte[i]) ^ multiply(e, a.byte[i]);
  }
  ret.byte[0] ^= de;
  return ret;
}
static inline aes_byte_share_t mulShareConst(aes_byte_share_t a, uint8_t b) {
  aes_byte_share_t ret;
  for (size_t i = 0; i < 64; i++) {
    ret.byte[i] = multiply(b, a.byte[i]);
  }
  return ret;
}

aes_byte_share_t mpc_aes_sbox(aes_byte_share_t s, randomTape_t* tapes, msgs_t* msgs) {

  // do the inversion
  // get random r, and a multiplication triple a,b,c, multiply r*s
  aes_byte_share_t r;
  uint8_t sr;
  do {
    r                  = getByteFromTapes(tapes);
    aes_byte_share_t a = getByteFromTapes(tapes);
    aes_byte_share_t b = getByteFromTapes(tapes);
    aes_byte_share_t c = getByteFromTapes(tapes);

	//uint8_t a1 = broadcastShare(a, msgs), b1 = broadcastShare(b, msgs), c1 = broadcastShare(c, msgs);
    //uint8_t r1 = broadcastShare(r, msgs), s1 = broadcastShare(s, msgs);
    //printf("%u,%u,%u,%u, %u,%u,%u, ", a1, b1, c1, multiply(a1,b1), r1, s1, multiply(r1, s1));
    //msgs->pos -= 5;
    // compute <s*r>
    aes_byte_share_t d_shares  = subShare(s, a);
    aes_byte_share_t e_shares  = subShare(r, b);
    uint8_t d                  = broadcastShare(d_shares, msgs);
    uint8_t e                  = broadcastShare(e_shares, msgs);
    //printf("%u,%u, ", d, e);
    aes_byte_share_t sr_shares = computeMul(a, b, c, d, e);
    sr                         = broadcastShare(sr_shares, msgs);
  } while (sr == 0);

  // invert sr locally
  uint8_t c3     = multiply(square(sr), sr);
  uint8_t c7     = multiply(square(c3), sr);
  uint8_t c63    = multiply(square(square(square(c7))), c7);
  uint8_t c127   = multiply(square(c63), sr);
  uint8_t sr_inv = square(c127);
  //printf("%u,%u,%u,%u,%u,%u, ", sr, c3, c7, c63, c127, sr_inv);

  // compute shares of sr_inv
  aes_byte_share_t sr_inv_shares = mulShareConst(r, sr_inv);

  uint8_t f[8];
  uint8_t h[8];
  aes_byte_share_t result = sr_inv_shares;
  int i;
  // finally, compute affine map per party, only doing the addition for the first party
  for (int j = 0; j < 64; j++) {
    uint8_t is_first = (j == 0);
    for (i = 0; i < 8; ++i)
      f[i] = 1 & (sr_inv_shares.byte[j] >> i);
    h[0]           = f[0] ^ f[4] ^ f[5] ^ f[6] ^ f[7] ^ is_first;
    h[1]           = f[1] ^ f[5] ^ f[6] ^ f[7] ^ f[0] ^ is_first;
    h[2]           = f[2] ^ f[6] ^ f[7] ^ f[0] ^ f[1];
    h[3]           = f[3] ^ f[7] ^ f[0] ^ f[1] ^ f[2];
    h[4]           = f[4] ^ f[0] ^ f[1] ^ f[2] ^ f[3];
    h[5]           = f[5] ^ f[1] ^ f[2] ^ f[3] ^ f[4] ^ is_first;
    h[6]           = f[6] ^ f[2] ^ f[3] ^ f[4] ^ f[5] ^ is_first;
    h[7]           = f[7] ^ f[3] ^ f[4] ^ f[5] ^ f[6];
    result.byte[j] = 0;
    for (i = 0; i < 8; ++i)
      result.byte[j] |= h[i] << i;
  }
  //uint8_t inp = broadcastShare(s, msgs);
  //msgs->pos--;
  //uint8_t res = broadcastShare(result, msgs);
  //msgs->pos--;
  //printf("%u->%u\n", inp, res);
  return result;
}

int mpc_aes_128(mzd_local_t* key, shares_aes_t* key_shares, randomTape_t* tapes, msgs_t* msgs,
                const uint8_t* plaintext, const uint8_t* pubKey, const picnic_instance_t* params) {

  (void)key;
  (void)params;
  aes_byte_share_t expanded[4][44];
  aes_byte_share_t state[4][4];
  aes_byte_share_t newstate[4][4];
  uint8_t ciphertext[16];
  uint8_t roundconstant;
  int i;
  int j;
  int r;
  int ret = 0;

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      expanded[i][j] = key_shares->shares[j * 4 + i];

  roundconstant = 1;
  for (j = 4; j < 44; ++j) {
    aes_byte_share_t temp[4];
    if (j % 4)
      for (i = 0; i < 4; ++i)
        temp[i] = expanded[i][j - 1];
    else {
      for (i = 0; i < 4; ++i) {
        temp[i] = mpc_aes_sbox(expanded[(i + 1) % 4][j - 1], tapes, msgs);
      }
      temp[0].byte[0] ^= roundconstant;
      roundconstant = xtime(roundconstant);
    }
    for (i = 0; i < 4; ++i)
      expanded[i][j] = subShare(temp[i], expanded[i][j - 4]);
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i) {
      state[i][j]         = expanded[i][j];
      state[i][j].byte[0] = plaintext[j * 4 + i] ^ expanded[i][j].byte[0];
    }

  for (r = 0; r < 10; ++r) {
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j) {
        newstate[i][j] = mpc_aes_sbox(state[i][j], tapes, msgs);
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] = newstate[i][(j + i) % 4];
    if (r < 9)
      for (j = 0; j < 4; ++j) {
        aes_byte_share_t a0 = state[0][j];
        aes_byte_share_t a1 = state[1][j];
        aes_byte_share_t a2 = state[2][j];
        aes_byte_share_t a3 = state[3][j];
        state[0][j] = subShare(subShare(subShare(mulShareConst(subShare(a0, a1), 2), a1), a2), a3);
        state[1][j] = subShare(subShare(subShare(mulShareConst(subShare(a1, a2), 2), a2), a3), a0);
        state[2][j] = subShare(subShare(subShare(mulShareConst(subShare(a2, a3), 2), a3), a0), a1);
        state[3][j] = subShare(subShare(subShare(mulShareConst(subShare(a3, a0), 2), a0), a1), a2);
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] = subShare(state[i][j], expanded[i][r * 4 + 4 + j]);
  }

  for (j = 0; j < 4; ++j)
    for (i = 0; i < 4; ++i)
      ciphertext[j * 4 + i] = broadcastShare(state[i][j], msgs);

  if (memcmp(ciphertext, pubKey, 16) != 0) {
#if 1 || !defined(NDEBUG)
    printf("%s: output does not match pubKey\n", __func__);
    printf("pubKey: ");
    print_hex(stdout, (uint8_t*)pubKey, 16);
    printf("\noutput: ");
    print_hex(stdout, (uint8_t*)ciphertext, 16);
    printf("\n");
#endif
    ret = -1;
  }

  return ret;
}
