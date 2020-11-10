#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../aes.h"

TEST_CASE("AES-128 KAT", "[aes]") {
  const std::vector<uint8_t> key = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                                    0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                                    0xff, 0xff, 0xff, 0xff};
  const std::vector<uint8_t> plaintext = {0x01, 0x01, 0x01, 0x01, 0x00, 0x00,
                                          0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                          0x00, 0x00, 0x00, 0x00};
  const std::vector<uint8_t> ciphertext_expected = {
      0x0b, 0x5a, 0x81, 0x4d, 0x95, 0x60, 0x1c, 0xc7,
      0xef, 0xe7, 0x12, 0x28, 0x3e, 0x05, 0xef, 0x8f};

  std::vector<uint8_t> ct;

  REQUIRE(AES128::aes_128(key, plaintext, ct) == true);
  REQUIRE(ct == ciphertext_expected);
}

TEST_CASE("AES-128 normal is equal to sbox-saving", "[aes]") {
  const std::vector<uint8_t> key = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                                    0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                                    0xff, 0xff, 0xff, 0xff};
  const std::vector<uint8_t> plaintext = {0x01, 0x01, 0x01, 0x01, 0x00, 0x00,
                                          0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                          0x00, 0x00, 0x00, 0x00};
  const std::vector<uint8_t> ciphertext_expected = {
      0x0b, 0x5a, 0x81, 0x4d, 0x95, 0x60, 0x1c, 0xc7,
      0xef, 0xe7, 0x12, 0x28, 0x3e, 0x05, 0xef, 0x8f};

  std::vector<uint8_t> ct;
  std::vector<uint8_t> ct2;

  REQUIRE(AES128::aes_128(key, plaintext, ct) == true);
  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> sbox_states =
      AES128::aes_128_with_sbox_output(key, plaintext, ct2);
  REQUIRE(ct == ciphertext_expected);
  REQUIRE(ct == ct2);
  REQUIRE(sbox_states.first.size() == AES128::NUM_SBOXES);
  REQUIRE(sbox_states.second.size() == AES128::NUM_SBOXES);
  // TODO hardcode values for this KAT
}