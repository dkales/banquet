#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../banquet.h"

TEST_CASE("Sign and verify a message", "[banquet]") {
  const char *message = "TestMessage";
  const banquet_instance_t &instance = banquet_instance_get(Banquet_L1_Param1);
  // banquet_keypair_t keypair = banquet_keygen(instance);
  const aes_block_t key = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                           0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
  const aes_block_t plaintext = {0x01, 0x01, 0x01, 0x01, 0x00, 0x00,
                                 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                 0x00, 0x00, 0x00, 0x00};
  const aes_block_t ciphertext_expected = {0x0b, 0x5a, 0x81, 0x4d, 0x95, 0x60,
                                           0x1c, 0xc7, 0xef, 0xe7, 0x12, 0x28,
                                           0x3e, 0x05, 0xef, 0x8f};

  banquet_keypair_t keypair;
  keypair.first = key;
  memcpy(keypair.second.data(), plaintext.data(), plaintext.size());
  memcpy(keypair.second.data() + plaintext.size(), ciphertext_expected.data(),
         ciphertext_expected.size());

  banquet_signature_t signature = banquet_sign(
      instance, keypair, (const uint8_t *)message, strlen(message));
  REQUIRE(signature.proofs.size() == instance.num_rounds);
  std::vector<uint8_t> serialized_signature =
      banquet_serialize_signature(instance, signature);
  std::cout << "signature length: " << serialized_signature.size()
            << " bytes\n";
}