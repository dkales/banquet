#pragma once

#include "gsl-lite.hpp"
#include "types.h"
#include <array>
#include <cstdint>
#include <vector>

namespace AES128 {
constexpr size_t NUM_SBOXES = 200;
constexpr size_t BLOCK_SIZE = 16;
constexpr size_t KEY_SIZE = 16;
constexpr size_t NUM_BLOCKS = 1;

bool aes_128(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintext_in,
             std::vector<uint8_t> &ciphertext_out);

std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
aes_128_with_sbox_output(const std::vector<uint8_t> &key_in,
                         const std::vector<uint8_t> &plaintext_in,
                         std::vector<uint8_t> &ciphertext_out);

void aes_128_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
                      const std::vector<gsl::span<uint8_t>> &t_shares,
                      const std::vector<uint8_t> &plaintext_in,
                      std::vector<gsl::span<uint8_t>> &ciphertext_shares_out,
                      std::vector<gsl::span<uint8_t>> &s_shares_out);

} // namespace AES128

namespace AES192 {
constexpr size_t NUM_SBOXES = 416;
constexpr size_t BLOCK_SIZE = 16;
constexpr size_t KEY_SIZE = 24;
constexpr size_t NUM_BLOCKS = 2;

bool aes_192(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintexts_in,
             std::vector<uint8_t> &ciphertexts_out);

std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
aes_192_with_sbox_output(const std::vector<uint8_t> &key_in,
                         const std::vector<uint8_t> &plaintext_in,
                         std::vector<uint8_t> &ciphertext_out);

void aes_192_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
                      const std::vector<gsl::span<uint8_t>> &t_shares,
                      const std::vector<uint8_t> &plaintext_in,
                      std::vector<gsl::span<uint8_t>> &ciphertext_shares_out,
                      std::vector<gsl::span<uint8_t>> &s_shares_out);

} // namespace AES192

namespace AES256 {
constexpr size_t NUM_SBOXES = 500;
constexpr size_t BLOCK_SIZE = 16;
constexpr size_t KEY_SIZE = 32;
constexpr size_t NUM_BLOCKS = 2;

bool aes_256(const std::vector<uint8_t> &key_in,
             const std::vector<uint8_t> &plaintexts_in,
             std::vector<uint8_t> &ciphertexts_out);

std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
aes_256_with_sbox_output(const std::vector<uint8_t> &key_in,
                         const std::vector<uint8_t> &plaintext_in,
                         std::vector<uint8_t> &ciphertext_out);

void aes_256_s_shares(const std::vector<gsl::span<uint8_t>> &key_in,
                      const std::vector<gsl::span<uint8_t>> &t_shares,
                      const std::vector<uint8_t> &plaintext_in,
                      std::vector<gsl::span<uint8_t>> &ciphertext_shares_out,
                      std::vector<gsl::span<uint8_t>> &s_shares_out);

} // namespace AES256