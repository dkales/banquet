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

static digest_t commit_to_party_seed(const seed_t &seed,
                                     const banquet_salt_t &salt, size_t rep_idx,
                                     size_t party_idx) {
  hash_context ctx;
  hash_init(&ctx, DIGEST_SIZE);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_idx);
  hash_update_uint16_le(&ctx, (uint16_t)party_idx);
  hash_update(&ctx, seed.data(), seed.size());
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}
static digest_t
phase_1_commitment(const banquet_instance_t &instance,
                   const banquet_salt_t &salt,
                   const std::vector<std::vector<digest_t>> &commitments,
                   const std::vector<aes_block_t> &key_deltas,
                   const std::vector<std::vector<uint8_t>> &t_deltas) {

  hash_context ctx;
  hash_init(&ctx, DIGEST_SIZE);
  hash_update(&ctx, salt.data(), salt.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      hash_update(&ctx, commitments[repetition][party].data(),
                  commitments[repetition][party].size());
      hash_update(&ctx, key_deltas[repetition].data(),
                  key_deltas[repetition].size());
      hash_update(&ctx, t_deltas[repetition].data(),
                  t_deltas[repetition].size());
    }
  }
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
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
  std::vector<std::vector<digest_t>> party_seed_commitments;

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    // generate seed tree for the N parties
    seed_t &repetition_seed = master_seeds[repetition];
    seed_trees.emplace_back(repetition_seed, instance.num_MPC_parties, salt,
                            repetition);

    // commit to each party's seed;
    std::vector<digest_t> current_party_seed_commitments;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      current_party_seed_commitments.push_back(
          commit_to_party_seed(seed_trees[repetition].get_leaf(party).value(),
                               salt, repetition, party));
    }
    party_seed_commitments.push_back(current_party_seed_commitments);

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
  std::vector<std::vector<aes_block_t>> rep_shared_keys;
  std::vector<aes_block_t> rep_key_deltas;
  std::vector<std::vector<std::vector<uint8_t>>> rep_shared_ts;
  std::vector<std::vector<uint8_t>> rep_t_deltas;

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

    // fix first share
    std::transform(std::begin(key_delta), std::end(key_delta),
                   std::begin(shared_key[0]), std::begin(shared_key[0]),
                   std::bit_xor<uint8_t>());

    rep_shared_keys.push_back(shared_key);
    rep_key_deltas.push_back(key_delta);
    // generate sharing of t values
    std::vector<std::vector<uint8_t>> shared_ts(instance.num_MPC_parties);
    std::vector<uint8_t> t_deltas(NUM_SBOXES_AES_128); // todo: copy initial ts
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      shared_ts[party].resize(NUM_SBOXES_AES_128);
      std::transform(std::begin(shared_ts[party]), std::end(shared_ts[party]),
                     std::begin(t_deltas), std::begin(t_deltas),
                     std::bit_xor<uint8_t>());
    }
    // fix first share
    std::transform(std::begin(t_deltas), std::end(t_deltas),
                   std::begin(shared_ts[0]), std::begin(shared_ts[0]),
                   std::bit_xor<uint8_t>());
    rep_shared_ts.push_back(shared_ts);
    rep_t_deltas.push_back(t_deltas);
  }
  // commit to salt, (all commitments of parties seeds, key_delta, t_delta) for
  // all repetitions
  digest_t sigma_1 = phase_1_commitment(instance, salt, party_seed_commitments,
                                        rep_key_deltas, rep_t_deltas);

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