#pragma once

#include "types.h"
#include <array>
#include <cstdint>
#include <vector>

constexpr size_t NUM_SBOXES_AES_128 = 200;

bool aes_128(const aes_block_t &key_in, const aes_block_t &plaintext_in,
             aes_block_t &ciphertext_out);

std::vector<std::pair<uint8_t, uint8_t>>
aes_128_with_sbox_output(const aes_block_t &key_in,
                         const aes_block_t &plaintext_in,
                         aes_block_t &ciphertext_out);

std::vector<std::vector<uint8_t>>
aes_128_s_shares(const std::vector<aes_block_t> &key_in,
                 const std::vector<std::vector<uint8_t>> &t_shares,
                 const aes_block_t &plaintext_in, aes_block_t &ciphertext_out);
