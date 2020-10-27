#include "aes.h"

int aes_128(const aes_block_t key_in, const aes_block_t plaintext_in, aes_block_t &ciphertext_out)
{
  return aes_128_old(key_in.data(), plaintext_in.data(), ciphertext_out.data());
}

static unsigned char multiply(unsigned int a, unsigned int b)
{
  unsigned char result = 0;

  for (int i = 0; i < 8; ++i)
  {
    uint8_t mask = -(b & 1);
    result ^= (a & mask);
    uint16_t mask2 = -((uint16_t)(a >> 7) & 1);
    a <<= 1;
    a ^= (0x11b & mask2);
    b >>= 1;
  }
  return result;
}

static unsigned char multiply2(unsigned int c, unsigned int d)
{
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

  for (i = 6; i >= 0; --i)
  {
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

static unsigned char square(unsigned char c)
{
  return multiply(c, c);
}

static unsigned char xtime(unsigned char c)
{
  return multiply(c, 2);
}

static unsigned char bytesub(unsigned char c)
{
  unsigned char c3 = multiply(square(c), c);
  unsigned char c7 = multiply(square(c3), c);
  unsigned char c63 = multiply(square(square(square(c7))), c7);
  unsigned char c127 = multiply(square(c63), c);
  unsigned char c254 = square(c127);
  unsigned char f[8];
  unsigned char h[8];
  unsigned char result = c254;
  int i;

  for (i = 0; i < 8; ++i)
    f[i] = 1 & (c254 >> i);
  h[0] = f[0] ^ f[4] ^ f[5] ^ f[6] ^ f[7] ^ 1;
  h[1] = f[1] ^ f[5] ^ f[6] ^ f[7] ^ f[0] ^ 1;
  h[2] = f[2] ^ f[6] ^ f[7] ^ f[0] ^ f[1];
  h[3] = f[3] ^ f[7] ^ f[0] ^ f[1] ^ f[2];
  h[4] = f[4] ^ f[0] ^ f[1] ^ f[2] ^ f[3];
  h[5] = f[5] ^ f[1] ^ f[2] ^ f[3] ^ f[4] ^ 1;
  h[6] = f[6] ^ f[2] ^ f[3] ^ f[4] ^ f[5] ^ 1;
  h[7] = f[7] ^ f[3] ^ f[4] ^ f[5] ^ f[6];
  result = 0;
  for (i = 0; i < 8; ++i)
    result |= h[i] << i;
  //printf("%u->%u\n", c, result);
  return result;
}

int aes_128_old(const uint8_t *key, const uint8_t *plaintext, uint8_t *ciphertext)
{
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
  for (j = 4; j < 44; ++j)
  {
    unsigned char temp[4];
    if (j % 4)
      for (i = 0; i < 4; ++i)
        temp[i] = expanded[i][j - 1];
    else
    {
      for (i = 0; i < 4; ++i)
      {
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

  for (r = 0; r < 10; ++r)
  {
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
      {
        if (state[i][j] == 0)
          return 1;
        newstate[i][j] = bytesub(state[i][j]);
      }
    for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
        state[i][j] = newstate[i][(j + i) % 4];
    if (r < 9)
      for (j = 0; j < 4; ++j)
      {
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

  return 0;
}
