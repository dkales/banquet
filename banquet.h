#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>

#include "banquet_instances.h"

/** Parameter set names */
enum banquet_params_t {
  PARAMETER_SET_INVALID = 0,
  Banquet_L1_Param1 = 1,
  PARAMETER_SET_MAX_INDEX = 2
};

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
banquet_keypair_t generate_banquet_key(banquet_instance_t instance);