#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include <NTL/GF2E.h>
using namespace NTL;

/* Prefix values for domain separation */
constexpr uint8_t HASH_PREFIX_0 = 0;
constexpr uint8_t HASH_PREFIX_1 = 1;
constexpr uint8_t HASH_PREFIX_2 = 2;
constexpr uint8_t HASH_PREFIX_3 = 3;
constexpr uint8_t HASH_PREFIX_4 = 4;
constexpr uint8_t HASH_PREFIX_5 = 5;

constexpr size_t SALT_SIZE = 32;
typedef std::array<uint8_t, SALT_SIZE> banquet_salt_t;

constexpr size_t SEED_SIZE = 16;
typedef std::array<uint8_t, SEED_SIZE> seed_t;

constexpr size_t DIGEST_SIZE = 32;
typedef std::array<uint8_t, DIGEST_SIZE> digest_t;

typedef std::pair<std::vector<seed_t>, size_t> reveal_list_t;

constexpr size_t BANQUET_PUBLICKEY_SIZE = 32;
constexpr size_t BANQUET_PRIVATEKEY_SIZE = 16;

typedef std::array<uint8_t, 16> aes_block_t;

/** Public key */
typedef std::array<uint8_t, BANQUET_PUBLICKEY_SIZE> banquet_publickey_t;

/** Private key */
typedef std::array<uint8_t, BANQUET_PRIVATEKEY_SIZE> banquet_privatekey_t;

typedef std::pair<banquet_privatekey_t, banquet_publickey_t> banquet_keypair_t;

struct banquet_repetition_proof_t {
  reveal_list_t reveallist;
  digest_t C_e;
  aes_block_t sk_delta;
  std::vector<uint8_t> t_delta;
  aes_block_t output_broadcast;
  std::vector<GF2E> P_delta;
  GF2E P_at_R;
  std::vector<GF2E> S_j_at_R;
  std::vector<GF2E> T_j_at_R;
};

struct banquet_signature_t {
  banquet_salt_t salt;
  digest_t h_1;
  digest_t h_2;
  digest_t h_3;
  std::vector<banquet_repetition_proof_t> proofs;
};