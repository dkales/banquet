#pragma once

#include "types.h"
#include <cstdint>
#include <array>

typedef std::array<uint8_t, 16> aes_block_t;

int aes_128(const aes_block_t &key_in, const aes_block_t &plaintext_in, aes_block_t &ciphertext_out);
int aes_128_old(const uint8_t *key, const uint8_t *plaintext, uint8_t *ciphertext);
