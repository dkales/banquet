#pragma once

#include "types.h"
#include <array>
#include <cstdint>
#include <vector>

typedef std::array<uint8_t, 16> aes_block_t;
constexpr size_t NUM_SBOXES_AES_128 = 200;

bool aes_128(const aes_block_t &key_in, const aes_block_t &plaintext_in,
             aes_block_t &ciphertext_out);

std::vector<std::pair<uint8_t, uint8_t>>
aes_128_with_sbox_output(const aes_block_t &key_in,
                         const aes_block_t &plaintext_in,
                         aes_block_t &ciphertext_out);
