#include "banquet.h"

#include "aes.h"
#include "kdf_shake.h"
#include "tape.h"
#include "tree.h"
#include "utils.h"
#include <algorithm>
#include <cassert>
#include <cstring>

#include <NTL/GF2EX.h>

banquet_keypair_t banquet_keygen(const banquet_instance_t &instance) {}

static void hash_update_GF2E(hash_context *ctx, const GF2E &element) {
  const GF2X &poly_rep = rep(element);
  std::vector<uint8_t> buffer(NumBytes(poly_rep));
  BytesFromGF2X(buffer.data(), poly_rep, buffer.size());

  hash_update(ctx, buffer.data(), buffer.size());
}

static std::pair<banquet_salt_t, std::vector<seed_t>>
generate_salt_and_seeds(const banquet_instance_t &instance,
                        const banquet_keypair_t &keypair,
                        const uint8_t *message, size_t message_len) {
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
  hash_init_prefix(&ctx, DIGEST_SIZE, HASH_PREFIX_1);
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

static std::vector<std::vector<GF2E>>
phase_1_expand(const banquet_instance_t &instance, const digest_t &h_1) {
  hash_context ctx;
  hash_init(&ctx, DIGEST_SIZE);
  hash_update(&ctx, h_1.data(), h_1.size());
  hash_final(&ctx);

  std::vector<std::vector<GF2E>> r_ejs;
  r_ejs.reserve(instance.num_rounds);
  for (size_t e = 0; e < instance.num_rounds; e++) {
    std::vector<GF2E> r_js;
    r_js.reserve(instance.m1);
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> r(instance.lambda);
      hash_squeeze(&ctx, r.data(), r.size());
      r_js.push_back(utils::GF2E_from_bytes(r));
    }
    r_ejs.push_back(r_js);
  }
  return r_ejs;
}
static digest_t
phase_2_commitment(const banquet_instance_t &instance,
                   const banquet_salt_t &salt, const digest_t &h_1,
                   const std::vector<std::vector<GF2E>> &P_deltas) {

  hash_context ctx;
  hash_init_prefix(&ctx, DIGEST_SIZE, HASH_PREFIX_2);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, h_1.data(), h_1.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t k = 0; k < instance.m2; k++) {
      hash_update_GF2E(&ctx, P_deltas[repetition][k]);
    }
  }
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

static std::vector<GF2E> phase_2_expand(const banquet_instance_t &instance,
                                        const digest_t &h_2,
                                        const vec_GF2E &forbidden_values) {
  hash_context ctx;
  hash_init(&ctx, DIGEST_SIZE);
  hash_update(&ctx, h_2.data(), h_2.size());
  hash_final(&ctx);

  std::vector<GF2E> R_es;
  for (size_t e = 0; e < instance.num_rounds; e++) {
    std::vector<uint8_t> R(instance.lambda);
    while (true) {
      hash_squeeze(&ctx, R.data(), R.size());
      //  check that R is not in {0,...m2-1}
      GF2E candidate_R = utils::GF2E_from_bytes(R);
      bool good = true;
      for (size_t k = 0; k < instance.m2; k++) {
        if (candidate_R == forbidden_values[k]) {
          good = false;
          break;
        }
      }
      if (good) {
        R_es.push_back(candidate_R);
        break;
      }
    }
  }
  return R_es;
}

static digest_t phase_3_commitment(
    const banquet_instance_t &instance, const banquet_salt_t &salt,
    const digest_t &h_2, const std::vector<GF2E> &c,
    const std::vector<std::vector<GF2E>> &c_shares,
    const std::vector<std::vector<GF2E>> &a,
    const std::vector<std::vector<std::vector<GF2E>>> &a_shares,
    const std::vector<std::vector<GF2E>> &b,
    const std::vector<std::vector<std::vector<GF2E>>> &b_shares) {

  hash_context ctx;
  hash_init_prefix(&ctx, DIGEST_SIZE, HASH_PREFIX_3);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, h_2.data(), h_2.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    hash_update_GF2E(&ctx, c[repetition]);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      hash_update_GF2E(&ctx, c_shares[repetition][party]);
    }
    for (size_t j = 0; j < instance.m1; j++) {
      hash_update_GF2E(&ctx, a[repetition][j]);
      hash_update_GF2E(&ctx, b[repetition][j]);
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        hash_update_GF2E(&ctx, a_shares[repetition][party][j]);
        hash_update_GF2E(&ctx, b_shares[repetition][party][j]);
      }
    }
  }
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

static std::vector<uint8_t> phase_3_expand(const banquet_instance_t &instance,
                                           const digest_t &h_3) {
  assert(instance.num_MPC_parties <= 256);
  // TODO assert(is_power_of_2(instance.num_MPC_parties));
  hash_context ctx;
  hash_init(&ctx, DIGEST_SIZE);
  hash_update(&ctx, h_3.data(), h_3.size());
  hash_final(&ctx);

  std::vector<uint8_t> opened_parties;
  uint8_t mask = (instance.num_MPC_parties - 1);
  for (size_t e = 0; e < instance.num_rounds; e++) {
    uint8_t party;
    hash_squeeze(&ctx, &party, sizeof(uint8_t));
    party = party & mask;
    opened_parties.push_back(party);
  }
  return opened_parties;
}

banquet_signature_t banquet_sign(const banquet_instance_t &instance,
                                 const banquet_keypair_t &keypair,
                                 const uint8_t *message, size_t message_len) {
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
  std::vector<std::vector<std::vector<uint8_t>>> rep_shared_s;
  std::vector<std::vector<std::vector<uint8_t>>> rep_shared_t;
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
    rep_shared_t.push_back(shared_ts);
    rep_t_deltas.push_back(t_deltas);
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 2: challenge the multiplications
  /////////////////////////////////////////////////////////////////////////////

  // commit to salt, (all commitments of parties seeds, key_delta, t_delta)
  // for all repetitions
  digest_t h_1 = phase_1_commitment(instance, salt, party_seed_commitments,
                                    rep_key_deltas, rep_t_deltas);

  // expand challenge hash to M * m1 values
  std::vector<std::vector<GF2E>> r_ejs = phase_1_expand(instance, h_1);

  /////////////////////////////////////////////////////////////////////////////
  // phase 3: commit to the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  // init modulus of extension field F_{2^{8\lambda}}
  utils::init_extension_field(instance);

  // a vector of the first m2+1 field elements for interpolation
  vec_GF2E x_values_for_interpolation_zero_to_m2 =
      utils::get_first_n_field_elements(instance.m2 + 1);
  vec_GF2E x_values_for_interpolation_zero_to_2m2 =
      utils::get_first_n_field_elements(2 * instance.m2 + 1);

  std::vector<std::vector<std::vector<GF2EX>>> S_eji(instance.num_rounds);
  std::vector<std::vector<std::vector<GF2EX>>> T_eji(instance.num_rounds);
  std::vector<GF2EX> P_e(instance.num_rounds);
  std::vector<std::vector<GF2EX>> P_ei(instance.num_rounds);
  std::vector<std::vector<GF2E>> P_deltas(instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    S_eji[repetition].resize(instance.num_MPC_parties);
    T_eji[repetition].resize(instance.num_MPC_parties);
    P_ei[repetition].resize(instance.num_MPC_parties);
    P_deltas[repetition].resize(instance.m2 + 1);

    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      S_eji[repetition][party].resize(instance.m1);
      T_eji[repetition][party].resize(instance.m1);
      // lift shares from F_{2^8} to F_{2^{8\lambda}}
      std::vector<GF2E> lifted_s;
      std::vector<GF2E> lifted_t;
      lifted_s.reserve(NUM_SBOXES_AES_128);
      lifted_t.reserve(NUM_SBOXES_AES_128);
      for (size_t idx = 0; idx < NUM_SBOXES_AES_128; idx++) {
        lifted_s.push_back(
            utils::lift_uint8_t(rep_shared_s[repetition][party][idx]));
        lifted_t.push_back(
            utils::lift_uint8_t(rep_shared_t[repetition][party][idx]));
      }
      std::vector<vec_GF2E> s_bar(instance.m1);
      std::vector<vec_GF2E> t_bar(instance.m1);

      for (size_t j = 0; j < instance.m1; j++) {
        // rearrange shares
        s_bar[j].SetLength(instance.m2 + 1);
        t_bar[j].SetLength(instance.m2 + 1);
        for (size_t k = 0; k < instance.m2; k++) {
          s_bar[j][k] = r_ejs[repetition][j] * lifted_s[j + instance.m1 * k];
          t_bar[j][k] = lifted_t[j + instance.m1 * k];
        }

        // sample additional random points
        std::vector<uint8_t> s_ej_bar(instance.lambda);
        random_tapes[repetition][party].squeeze_bytes(s_ej_bar.data(),
                                                      s_ej_bar.size());

        std::vector<uint8_t> t_ej_bar(instance.lambda);
        random_tapes[repetition][party].squeeze_bytes(t_ej_bar.data(),
                                                      t_ej_bar.size());
        s_bar[j][instance.m2] = utils::GF2E_from_bytes(s_ej_bar);
        t_bar[j][instance.m2] = utils::GF2E_from_bytes(t_ej_bar);

        // interpolate polynomials S_ej^i and T_ej^i
        S_eji[repetition][party][j] =
            interpolate(x_values_for_interpolation_zero_to_m2, s_bar[j]);
        T_eji[repetition][party][j] =
            interpolate(x_values_for_interpolation_zero_to_m2, t_bar[j]);
      }
    }

    // compute product polynomial P_e
    GF2EX P;
    for (size_t j = 0; j < instance.m1; j++) {
      GF2EX S_sum;
      GF2EX T_sum;

      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        S_sum += S_eji[repetition][party][j];
        T_sum += T_eji[repetition][party][j];
      }

      P += S_sum * T_sum;
    }
    P_e[repetition] = P;

    // compute sharing of P
    std::vector<vec_GF2E> P_shares(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // first m2 points: first party = sum of r_e,j, other parties = 0
      P_shares[party].SetLength(2 * instance.m2 + 1);
      if (party == 0) {
        GF2E sum_r;
        for (size_t j = 0; j < instance.m1; j++) {
          sum_r += r_ejs[repetition][j];
        }
        for (size_t k = 0; k < instance.m2; k++) {
          P_shares[party][k] = sum_r;
        }
      } else {
        for (size_t k = 0; k < instance.m2; k++) {
          P_shares[party][k] = GF2E::zero();
        }
      }

      // second m2+1 points: sample from random tape
      for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
        std::vector<uint8_t> P_k_share(instance.lambda);
        random_tapes[repetition][party].squeeze_bytes(P_k_share.data(),
                                                      P_k_share.size());
        P_shares[party][k] = utils::GF2E_from_bytes(P_k_share);
      }
    }
    for (size_t k = instance.m2; k < 2 * instance.m2; k++) {
      // calculate offset
      GF2E k_element = x_values_for_interpolation_zero_to_2m2[k];
      GF2E P_at_k_delta = eval(P, k_element);
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        P_at_k_delta -= P_shares[party][k];
      }
      P_deltas[repetition][k] = P_at_k_delta;
      // adjust first share
      P_shares[0][k] += P_at_k_delta;
    }
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // iterpolate polynomial P_e^1 from 2m+1 points
      P_ei[repetition][party] =
          interpolate(x_values_for_interpolation_zero_to_2m2, P_shares[party]);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 4: challenge the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  digest_t h_2 = phase_2_commitment(instance, salt, h_1, P_deltas);

  // expand challenge hash to M values

  vec_GF2E forbidden_challenge_values =
      utils::get_first_n_field_elements(instance.m2);
  std::vector<GF2E> R_es =
      phase_2_expand(instance, h_2, forbidden_challenge_values);

  /////////////////////////////////////////////////////////////////////////////
  // phase 5: commit to the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  std::vector<GF2E> c(instance.num_rounds);
  std::vector<std::vector<GF2E>> c_shares(instance.num_rounds);
  std::vector<std::vector<GF2E>> a(instance.num_rounds);
  std::vector<std::vector<std::vector<GF2E>>> a_shares(instance.num_rounds);
  std::vector<std::vector<GF2E>> b(instance.num_rounds);
  std::vector<std::vector<std::vector<GF2E>>> b_shares(instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    c_shares[repetition].resize(instance.num_MPC_parties);
    a[repetition].resize(instance.m1);
    b[repetition].resize(instance.m1);
    a_shares[repetition].resize(instance.num_MPC_parties);
    b_shares[repetition].resize(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      a_shares[repetition][party].resize(instance.m1);
      b_shares[repetition][party].resize(instance.m1);
      for (size_t j = 0; j < instance.m1; j++) {
        // compute a_ej^i and b_ej^i
        a_shares[repetition][party][j] =
            eval(S_eji[repetition][party][j], R_es[repetition]);
        b_shares[repetition][party][j] =
            eval(T_eji[repetition][party][j], R_es[repetition]);
      }
      // compute c_e^i
      c_shares[repetition][party] =
          eval(P_ei[repetition][party], R_es[repetition]);
    }
    // open c_e and a,b values
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      c[repetition] += c_shares[repetition][party];
      for (size_t j = 0; j < instance.m1; j++) {
        a[repetition][j] += a_shares[repetition][party][j];
        b[repetition][j] += b_shares[repetition][party][j];
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 6: challenge the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  digest_t h_3 = phase_3_commitment(instance, salt, h_2, c, c_shares, a,
                                    a_shares, b, b_shares);

  std::vector<uint8_t> missing_parties = phase_3_expand(instance, h_3);

  /////////////////////////////////////////////////////////////////////////////
  // phase 7: Open the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////
  std::vector<reveal_list_t> seeds;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    seeds.push_back(
        seed_trees[repetition].reveal_all_but(missing_parties[repetition]));
  }
  // build signature
  std::vector<banquet_repetition_proof_t> proofs;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    size_t missing_party = missing_parties[repetition];
    banquet_repetition_proof_t proof{
        seeds[repetition],
        party_seed_commitments[repetition][missing_party],
        rep_key_deltas[repetition],
        rep_t_deltas[repetition],
        P_deltas[repetition],
        c[repetition],
        a[repetition],
        b[repetition],
    };
    proofs.push_back(proof);
  }

  banquet_signature_t signature{salt, h_1, h_2, h_3, proofs};

  return signature;
}