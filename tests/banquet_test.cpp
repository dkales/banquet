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
TEST_CASE("Serialization and Deserialization", "[banquet]") {
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
  banquet_signature_t signature2 =
      banquet_deserialize_signature(instance, serialized_signature);
  REQUIRE(signature.salt == signature2.salt);
  REQUIRE(signature.h_1 == signature2.h_1);
  REQUIRE(signature.h_2 == signature2.h_2);
  REQUIRE(signature.h_3 == signature2.h_3);
  for (size_t i = 0; i < instance.num_rounds; i++) {
    banquet_repetition_proof_t &proof1 = signature.proofs[i];
    banquet_repetition_proof_t &proof2 = signature2.proofs[i];
    REQUIRE(proof1.reveallist == proof2.reveallist);
    REQUIRE(proof1.C_e == proof2.C_e);
    REQUIRE(proof1.sk_delta == proof2.sk_delta);
    REQUIRE(proof1.t_delta == proof2.t_delta);
    REQUIRE(proof1.P_delta == proof2.P_delta);
    REQUIRE(proof1.P_at_R == proof2.P_at_R);
    REQUIRE(proof1.S_j_at_R == proof2.S_j_at_R);
    REQUIRE(proof1.T_j_at_R == proof2.T_j_at_R);
  }
}