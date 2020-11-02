#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../banquet.h"

TEST_CASE("Sign and verify a message", "[banquet]") {
  const char *message = "TestMessage";
  const banquet_instance_t &instance = banquet_instance_get(Banquet_L1_Param1);
  banquet_keypair_t keypair = banquet_keygen(instance);
  banquet_signature_t signature = banquet_sign(
      instance, keypair, (const uint8_t *)message, strlen(message));
  REQUIRE(signature.proofs.size() == instance.num_rounds);
}