#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "banquet_instances.h"

/* Prefix values for domain separation */
constexpr uint8_t HASH_PREFIX_0 = 0;
constexpr uint8_t HASH_PREFIX_1 = 1;
constexpr uint8_t HASH_PREFIX_2 = 2;
constexpr uint8_t HASH_PREFIX_3 = 3;
constexpr uint8_t HASH_PREFIX_4 = 4;
constexpr uint8_t HASH_PREFIX_5 = 5;

constexpr size_t SALT_SIZE = 32;
typedef std::array<uint8_t, SALT_SIZE> banquet_salt_t;
constexpr size_t MAX_DIGEST_SIZE = 32;
constexpr size_t BANQUET_PUBLICKEY_SIZE = 32;
constexpr size_t BANQUET_PRIVATEKEY_SIZE = 16;

/** Public key */
typedef std::array<uint8_t, BANQUET_PUBLICKEY_SIZE> banquet_publickey_t;

/** Private key */
typedef std::array<uint8_t, BANQUET_PRIVATEKEY_SIZE> banquet_privatekey_t;

typedef std::pair<banquet_privatekey_t, banquet_publickey_t> banquet_keypair_t;

// crypto api
banquet_keypair_t banquet_keygen(const banquet_instance_t &instance);

std::vector<uint8_t> sign_message(const banquet_instance_t &instance,
                                  const banquet_keypair_t &keypair,
                                  uint8_t *message, size_t message_len);