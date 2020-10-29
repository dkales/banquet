#include "banquet.h"

#include "aes.h"
#include "kdf_shake.h"
#include "tape.h"
#include "tree.h"
#include <algorithm>
#include <cassert>
#include <cstring>

banquet_keypair_t banquet_keygen(const banquet_instance_t &instance) {}

static std::pair<banquet_salt_t, std::vector<seed_t>>
generate_salt_and_seeds(const banquet_instance_t &instance,
                        const banquet_keypair_t &keypair, uint8_t *message,
                        size_t message_len) {
  // salt, seed_1, ..., seed_r = H(instance||sk||pk||m)
  hash_context ctx;
  hash_init(&ctx, DIGEST_SIZE);
  hash_update_uint16_le(&ctx, (uint16_t)instance.params);
  hash_update(&ctx, keypair.first.data(), keypair.first.size());
  hash_update(&ctx, keypair.second.data(), keypair.second.size());
  hash_update(&ctx, message, message_len);
  hash_final(&ctx);

  banquet_salt_t salt;
  hash_squeeze(&ctx, salt.data(), salt.size());
  std::vector<seed_t> seeds(instance.num_rounds);
  for (seed_t &s : seeds) {
    hash_squeeze(&ctx, s.data(), s.size());
  }
  return std::make_pair(salt, seeds);
}

std::vector<uint8_t> banquet_sign(const banquet_instance_t &instance,
                                  const banquet_keypair_t &keypair,
                                  uint8_t *message, size_t message_len) {

  // grab aes key, pt and ct
  aes_block_t key = keypair.first;
  aes_block_t pt, ct, ct2;
  memcpy(pt.data(), keypair.second.data(), pt.size());
  memcpy(ct.data(), keypair.second.data() + pt.size(), ct.size());

  // get sbox inputs and outputs for aes evaluation
  std::vector<std::pair<uint8_t, uint8_t>> sbox_pairs =
      aes_128_with_sbox_output(key, pt, ct2);
  // sanity check, incoming keypair is valid
  assert(ct == ct2);

  // generate salt and master seeds for each repetition
  auto [salt, master_seeds] =
      generate_salt_and_seeds(instance, keypair, message, message_len);

  // do parallel repetitions
  // create seed trees and random tapes
  std::vector<SeedTree> seed_trees;
  std::vector<std::vector<RandomTape>> random_tapes;

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    // generate seed tree for the N parties
    seed_t &repetition_seed = master_seeds[repetition];
    seed_trees.emplace_back(repetition_seed, instance.num_MPC_parties, salt,
                            repetition);

    // create random tape for each party
    std::vector<RandomTape> party_tapes;
    party_tapes.reserve(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      party_tapes.emplace_back(seed_trees[repetition].get_leaf(party).value(),
                               salt, repetition, party);
    }
    random_tapes.push_back(party_tapes);
  }
  /////////////////////////////////////////////////////////////////////////////
  // phase 1: commit to executions of AES
  /////////////////////////////////////////////////////////////////////////////
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    // generate sharing of secret key
    std::vector<aes_block_t> shared_key(instance.num_MPC_parties);
    aes_block_t key_delta = key;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      random_tapes[repetition][party].squeeze_bytes(shared_key[party].data(),
                                                    shared_key[party].size());
      std::transform(std::begin(shared_key[party]), std::end(shared_key[party]),
                     std::begin(key_delta), std::begin(key_delta),
                     std::bit_xor<uint8_t>());
    }
    // generate sharing of t values
    std::vector<std::vector<uint8_t>> shared_ts(instance.num_MPC_parties);
    std::vector<uint8_t> t_deltas(NUM_SBOXES_AES_128);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      shared_ts[party].resize(NUM_SBOXES_AES_128);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 2: challenge the multiplications
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // phase 3: commit to the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // phase 4: challenge the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // phase 5: commit to the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // phase 6: challenge the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // phase 7: Open the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  // serialize signature
}