#include "banquet.h"

#include "aes.h"
#include "tape.h"
#include "tree.h"
#include "utils.h"
#include <algorithm>
#include <cassert>
#include <cstring>

extern "C" {
#include "kdf_shake.h"
#include "randomness.h"
}

#include <NTL/GF2EX.h>

namespace {
inline std::vector<uint8_t> GF2E_to_bytes(const banquet_instance_t &instance,
                                          const GF2E &element) {
  const GF2X &poly_rep = rep(element);
  std::vector<uint8_t> buffer(instance.lambda);
  BytesFromGF2X(buffer.data(), poly_rep, buffer.size());
  return buffer;
}

inline void hash_update_GF2E(hash_context *ctx,
                             const banquet_instance_t &instance,
                             const GF2E &element) {
  std::vector<uint8_t> buffer = GF2E_to_bytes(instance, element);
  hash_update(ctx, buffer.data(), buffer.size());
}

std::pair<banquet_salt_t, std::vector<seed_t>>
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

digest_t commit_to_party_seed(const seed_t &seed, const banquet_salt_t &salt,
                              size_t rep_idx, size_t party_idx) {
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

digest_t phase_1_commitment(
    const banquet_instance_t &instance, const banquet_salt_t &salt,
    const banquet_publickey_t &pk, const uint8_t *message, size_t message_len,
    const std::vector<std::vector<digest_t>> &commitments,
    const std::vector<aes_block_t> &key_deltas,
    const std::vector<std::vector<uint8_t>> &t_deltas,
    const std::vector<std::vector<aes_block_t>> &output_broadcasts) {

  hash_context ctx;
  hash_init_prefix(&ctx, DIGEST_SIZE, HASH_PREFIX_1);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, pk.data(), pk.size());
  hash_update(&ctx, message, message_len);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      hash_update(&ctx, commitments[repetition][party].data(),
                  commitments[repetition][party].size());
      hash_update(&ctx, output_broadcasts[repetition][party].data(),
                  output_broadcasts[repetition][party].size());
    }
    hash_update(&ctx, key_deltas[repetition].data(),
                key_deltas[repetition].size());
    hash_update(&ctx, t_deltas[repetition].data(), t_deltas[repetition].size());
  }
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<std::vector<GF2E>>
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

digest_t phase_2_commitment(const banquet_instance_t &instance,
                            const banquet_salt_t &salt, const digest_t &h_1,
                            const std::vector<std::vector<GF2E>> &P_deltas) {

  hash_context ctx;
  hash_init_prefix(&ctx, DIGEST_SIZE, HASH_PREFIX_2);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, h_1.data(), h_1.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t k = 0; k < instance.m2; k++) {
      hash_update_GF2E(&ctx, instance, P_deltas[repetition][k]);
    }
  }
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<GF2E> phase_2_expand(const banquet_instance_t &instance,
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

digest_t phase_3_commitment(
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
    hash_update_GF2E(&ctx, instance, c[repetition]);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      hash_update_GF2E(&ctx, instance, c_shares[repetition][party]);
    }
    for (size_t j = 0; j < instance.m1; j++) {
      hash_update_GF2E(&ctx, instance, a[repetition][j]);
      hash_update_GF2E(&ctx, instance, b[repetition][j]);
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        hash_update_GF2E(&ctx, instance, a_shares[repetition][party][j]);
        hash_update_GF2E(&ctx, instance, b_shares[repetition][party][j]);
      }
    }
  }
  hash_final(&ctx);

  digest_t commitment;
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<uint8_t> phase_3_expand(const banquet_instance_t &instance,
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
} // namespace

banquet_keypair_t banquet_keygen(const banquet_instance_t &instance) {
  aes_block_t key, pt, ct;

  while (true) {
    rand_bytes(key.data(), key.size());
    rand_bytes(pt.data(), pt.size());
    if (aes_128(key, pt, ct)) {
      break;
    }
  }
  banquet_keypair_t keypair;
  memcpy(keypair.first.data(), key.data(), key.size());
  memcpy(keypair.second.data(), pt.data(), pt.size());
  memcpy(keypair.second.data() + pt.size(), ct.data(), ct.size());
  return keypair;
}

banquet_signature_t banquet_sign(const banquet_instance_t &instance,
                                 const banquet_keypair_t &keypair,
                                 const uint8_t *message, size_t message_len) {
  // init modulus of extension field F_{2^{8\lambda}}
  utils::init_extension_field(instance);

  // grab aes key, pt and ct
  aes_block_t key = keypair.first;
  aes_block_t pt, ct, ct2;
  memcpy(pt.data(), keypair.second.data(), pt.size());
  memcpy(ct.data(), keypair.second.data() + pt.size(), ct.size());

  // get sbox inputs and outputs for aes evaluation
  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> sbox_pairs =
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
  std::vector<std::vector<aes_block_t>> rep_output_broadcasts;
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
    std::vector<std::vector<uint8_t>> shared_t(instance.num_MPC_parties);
    std::vector<uint8_t> t_deltas = sbox_pairs.second;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      shared_t[party].resize(NUM_SBOXES_AES_128);
      random_tapes[repetition][party].squeeze_bytes(shared_t[party].data(),
                                                    shared_t[party].size());
      std::transform(std::begin(shared_t[party]), std::end(shared_t[party]),
                     std::begin(t_deltas), std::begin(t_deltas),
                     std::bit_xor<uint8_t>());
    }
    // fix first share
    std::transform(std::begin(t_deltas), std::end(t_deltas),
                   std::begin(shared_t[0]), std::begin(shared_t[0]),
                   std::bit_xor<uint8_t>());
    rep_shared_t.push_back(shared_t);
    rep_t_deltas.push_back(t_deltas);

    // get shares of sbox inputs by executing MPC AES
    std::vector<aes_block_t> ct_shares;
    std::vector<std::vector<uint8_t>> shared_s =
        aes_128_s_shares(shared_key, shared_t, pt, ct_shares);

    // sanity check, mpc execution = plain one
    aes_block_t ct_check;
    memset(ct_check.data(), 0, ct_check.size());
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      std::transform(std::begin(ct_shares[party]), std::end(ct_shares[party]),
                     std::begin(ct_check), std::begin(ct_check),
                     std::bit_xor<uint8_t>());
    }

    assert(ct == ct_check);
    rep_shared_s.push_back(shared_s);
    rep_output_broadcasts.push_back(ct_shares);
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 2: challenge the multiplications
  /////////////////////////////////////////////////////////////////////////////

  // commit to salt, (all commitments of parties seeds, key_delta, t_delta)
  // for all repetitions
  digest_t h_1 =
      phase_1_commitment(instance, salt, keypair.second, message, message_len,
                         party_seed_commitments, rep_key_deltas, rep_t_deltas,
                         rep_output_broadcasts);

  // expand challenge hash to M * m1 values
  std::vector<std::vector<GF2E>> r_ejs = phase_1_expand(instance, h_1);

  /////////////////////////////////////////////////////////////////////////////
  // phase 3: commit to the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  // a vector of the first m2+1 field elements for interpolation
  vec_GF2E x_values_for_interpolation_zero_to_m2 =
      utils::get_first_n_field_elements(instance.m2 + 1);
  std::vector<GF2EX> precomputation_for_zero_to_m2 =
      utils::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_m2);
  vec_GF2E x_values_for_interpolation_zero_to_2m2 =
      utils::get_first_n_field_elements(2 * instance.m2 + 1);
  std::vector<GF2EX> precomputation_for_zero_to_2m2 =
      utils::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_2m2);

  std::vector<std::vector<std::vector<vec_GF2E>>> s_prime(instance.num_rounds);
  std::vector<std::vector<std::vector<vec_GF2E>>> t_prime(instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> S_eji(instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> T_eji(instance.num_rounds);
  std::vector<GF2EX> P_e(instance.num_rounds);
  std::vector<std::vector<GF2EX>> P_ei(instance.num_rounds);
  std::vector<std::vector<GF2E>> P_deltas(instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    s_prime[repetition].resize(instance.num_MPC_parties);
    t_prime[repetition].resize(instance.num_MPC_parties);
    // S_eji[repetition].resize(instance.num_MPC_parties);
    // T_eji[repetition].resize(instance.num_MPC_parties);
    P_ei[repetition].resize(instance.num_MPC_parties);
    P_deltas[repetition].resize(instance.m2 + 1);

    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      s_prime[repetition][party].resize(instance.m1);
      t_prime[repetition][party].resize(instance.m1);
      // S_eji[repetition][party].resize(instance.m1);
      // T_eji[repetition][party].resize(instance.m1);
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

      for (size_t j = 0; j < instance.m1; j++) {
        vec_GF2E s_bar;
        vec_GF2E t_bar;
        // rearrange shares
        s_bar.SetLength(instance.m2 + 1);
        t_bar.SetLength(instance.m2 + 1);
        for (size_t k = 0; k < instance.m2; k++) {
          s_bar[k] = r_ejs[repetition][j] * lifted_s[j + instance.m1 * k];
          t_bar[k] = lifted_t[j + instance.m1 * k];
        }

        // sample additional random points
        std::vector<uint8_t> s_ej_bar(instance.lambda);
        random_tapes[repetition][party].squeeze_bytes(s_ej_bar.data(),
                                                      s_ej_bar.size());

        std::vector<uint8_t> t_ej_bar(instance.lambda);
        random_tapes[repetition][party].squeeze_bytes(t_ej_bar.data(),
                                                      t_ej_bar.size());
        s_bar[instance.m2] = utils::GF2E_from_bytes(s_ej_bar);
        t_bar[instance.m2] = utils::GF2E_from_bytes(t_ej_bar);

        // interpolate polynomials S_ej^i and T_ej^i
        // S_eji[repetition][party][j] = utils::interpolate_with_precomputation(
        // precomputation_for_zero_to_m2, s_bar);
        // T_eji[repetition][party][j] =
        // utils::interpolate_with_precomputation(
        // precomputation_for_zero_to_m2, t_bar);
        s_prime[repetition][party][j] = s_bar;
        t_prime[repetition][party][j] = t_bar;
      }
    }

    // compute product polynomial P_e
    GF2EX P;
    for (size_t j = 0; j < instance.m1; j++) {
      vec_GF2E s_prime_sum;
      vec_GF2E t_prime_sum;
      s_prime_sum.SetLength(instance.m2 + 1);
      t_prime_sum.SetLength(instance.m2 + 1);
      GF2EX S_sum;
      GF2EX T_sum;

      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        // S_sum += S_eji[repetition][party][j];
        // T_sum += T_eji[repetition][party][j];
        s_prime_sum += s_prime[repetition][party][j];
        t_prime_sum += t_prime[repetition][party][j];
      }
      S_sum = utils::interpolate_with_precomputation(
          precomputation_for_zero_to_m2, s_prime_sum);
      T_sum = utils::interpolate_with_precomputation(
          precomputation_for_zero_to_m2, t_prime_sum);

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
    for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
      // calculate offset
      GF2E k_element = x_values_for_interpolation_zero_to_2m2[k];
      GF2E P_at_k_delta = eval(P, k_element);
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        P_at_k_delta -= P_shares[party][k];
      }
      P_deltas[repetition][k - instance.m2] = P_at_k_delta;
      // adjust first share
      P_shares[0][k] += P_at_k_delta;
    }
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // iterpolate polynomial P_e^1 from 2m+1 points
      P_ei[repetition][party] = utils::interpolate_with_precomputation(
          precomputation_for_zero_to_2m2, P_shares[party]);
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
    vec_GF2E lagrange_polys_evaluated_at_Re;
    lagrange_polys_evaluated_at_Re.SetLength(instance.m2 + 1);
    for (size_t k = 0; k < instance.m2 + 1; k++) {
      lagrange_polys_evaluated_at_Re[k] =
          eval(precomputation_for_zero_to_m2[k], R_es[repetition]);
    }

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
        // a_shares[repetition][party][j] =
        // eval(S_eji[repetition][party][j], R_es[repetition]);
        // b_shares[repetition][party][j] =
        // eval(T_eji[repetition][party][j], R_es[repetition]);
        a_shares[repetition][party][j] =
            lagrange_polys_evaluated_at_Re * s_prime[repetition][party][j];
        b_shares[repetition][party][j] =
            lagrange_polys_evaluated_at_Re * t_prime[repetition][party][j];
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
        rep_output_broadcasts[repetition][missing_party],
        P_deltas[repetition],
        c[repetition],
        a[repetition],
        b[repetition],
    };
    // sanity check c = sum_j a*b
    // GF2E accum;
    // clear(accum);
    // for (size_t j = 0; j < instance.m1; j++) {
    //   accum += a[repetition][j] * b[repetition][j];
    // }
    // if (accum != c[repetition])
    //   throw std::runtime_error("something wrong here");
    proofs.push_back(proof);
  }

  banquet_signature_t signature{salt, h_1, h_2, h_3, proofs};

  return signature;
}

bool banquet_verify(const banquet_instance_t &instance,
                    const banquet_publickey_t &pk,
                    const banquet_signature_t &signature,
                    const uint8_t *message, size_t message_len) {

  // init modulus of extension field F_{2^{8\lambda}}
  utils::init_extension_field(instance);

  aes_block_t pt, ct;
  memcpy(pt.data(), pk.data(), pt.size());
  memcpy(ct.data(), pk.data() + pt.size(), ct.size());

  // do parallel repetitions
  // create seed trees and random tapes
  std::vector<SeedTree> seed_trees;
  std::vector<std::vector<RandomTape>> random_tapes;
  std::vector<std::vector<digest_t>> party_seed_commitments;

  // recompute h_2
  std::vector<std::vector<GF2E>> P_deltas;
  for (const banquet_repetition_proof_t &proof : signature.proofs) {
    P_deltas.push_back(proof.P_delta);
  }
  digest_t h_2 =
      phase_2_commitment(instance, signature.salt, signature.h_1, P_deltas);

  // compute challenges based on hashes
  // h1 expansion
  std::vector<std::vector<GF2E>> r_ejs =
      phase_1_expand(instance, signature.h_1);
  // h2 expansion
  vec_GF2E forbidden_challenge_values =
      utils::get_first_n_field_elements(instance.m2);
  std::vector<GF2E> R_es =
      phase_2_expand(instance, h_2, forbidden_challenge_values);
  // h3 expansion already happened in deserialize to get missing parties
  std::vector<uint8_t> missing_parties =
      phase_3_expand(instance, signature.h_3);

  // rebuild SeedTrees
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];
    // regenerate generate seed tree for the N parties (except the missing
    // one)
    if (missing_parties[repetition] != proof.reveallist.second)
      throw std::runtime_error(
          "modified signature between deserialization and verify");
    seed_trees.emplace_back(proof.reveallist, instance.num_MPC_parties,
                            signature.salt, repetition);
    // commit to each party's seed, fill up missing one with data from proof
    std::vector<digest_t> current_party_seed_commitments;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        current_party_seed_commitments.push_back(
            commit_to_party_seed(seed_trees[repetition].get_leaf(party).value(),
                                 signature.salt, repetition, party));
      } else {
        current_party_seed_commitments.push_back(proof.C_e);
      }
    }
    party_seed_commitments.push_back(current_party_seed_commitments);

    // create random tape for each party, dummy one for missing party
    std::vector<RandomTape> party_tapes;
    party_tapes.reserve(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        party_tapes.emplace_back(seed_trees[repetition].get_leaf(party).value(),
                                 signature.salt, repetition, party);
      } else {
        party_tapes.emplace_back(seed_t{}, signature.salt, repetition, party);
      }
    }
    random_tapes.push_back(party_tapes);
  }
  /////////////////////////////////////////////////////////////////////////////
  // recompute commitments to executions of AES
  /////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<aes_block_t>> rep_shared_keys;
  std::vector<std::vector<std::vector<uint8_t>>> rep_shared_s;
  std::vector<std::vector<std::vector<uint8_t>>> rep_shared_t;
  std::vector<std::vector<aes_block_t>> rep_output_broadcasts;

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];

    // generate sharing of secret key
    std::vector<aes_block_t> shared_key(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      random_tapes[repetition][party].squeeze_bytes(shared_key[party].data(),
                                                    shared_key[party].size());
    }

    // fix first share
    std::transform(std::begin(proof.sk_delta), std::end(proof.sk_delta),
                   std::begin(shared_key[0]), std::begin(shared_key[0]),
                   std::bit_xor<uint8_t>());

    rep_shared_keys.push_back(shared_key);
    // generate sharing of t values
    std::vector<std::vector<uint8_t>> shared_t(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      shared_t[party].resize(NUM_SBOXES_AES_128);
      random_tapes[repetition][party].squeeze_bytes(shared_t[party].data(),
                                                    shared_t[party].size());
    }
    // fix first share
    std::transform(std::begin(proof.t_delta), std::end(proof.t_delta),
                   std::begin(shared_t[0]), std::begin(shared_t[0]),
                   std::bit_xor<uint8_t>());
    rep_shared_t.push_back(shared_t);

    // get shares of sbox inputs by executing MPC AES
    std::vector<aes_block_t> ct_shares;
    std::vector<std::vector<uint8_t>> shared_s =
        aes_128_s_shares(shared_key, shared_t, pt, ct_shares);

    // get missing output broadcast from proof
    ct_shares[missing_parties[repetition]] = proof.output_broadcast;
    // sanity check, mpc execution = plain one
    aes_block_t ct_check;
    memset(ct_check.data(), 0, ct_check.size());
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      std::transform(std::begin(ct_shares[party]), std::end(ct_shares[party]),
                     std::begin(ct_check), std::begin(ct_check),
                     std::bit_xor<uint8_t>());
    }
    // check MPC execution is correct
    if (memcmp(ct_check.data(), ct.data(), ct.size()) != 0) {
      return false;
    }

    rep_shared_s.push_back(shared_s);
    rep_output_broadcasts.push_back(ct_shares);
  }

  /////////////////////////////////////////////////////////////////////////////
  // recompute shares of polynomials
  /////////////////////////////////////////////////////////////////////////////
  // a vector of the first m2+1 field elements for interpolation
  vec_GF2E x_values_for_interpolation_zero_to_m2 =
      utils::get_first_n_field_elements(instance.m2 + 1);
  std::vector<GF2EX> precomputation_for_zero_to_m2 =
      utils::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_m2);
  vec_GF2E x_values_for_interpolation_zero_to_2m2 =
      utils::get_first_n_field_elements(2 * instance.m2 + 1);
  std::vector<GF2EX> precomputation_for_zero_to_2m2 =
      utils::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_2m2);

  std::vector<std::vector<std::vector<vec_GF2E>>> s_prime(instance.num_rounds);
  std::vector<std::vector<std::vector<vec_GF2E>>> t_prime(instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> S_eji(instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> T_eji(instance.num_rounds);
  std::vector<GF2EX> P_e(instance.num_rounds);
  std::vector<std::vector<GF2EX>> P_ei(instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];
    // S_eji[repetition].resize(instance.num_MPC_parties);
    // T_eji[repetition].resize(instance.num_MPC_parties);
    s_prime[repetition].resize(instance.num_MPC_parties);
    t_prime[repetition].resize(instance.num_MPC_parties);
    P_ei[repetition].resize(instance.num_MPC_parties);
    P_deltas[repetition].resize(instance.m2 + 1);

    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        // S_eji[repetition][party].resize(instance.m1);
        // T_eji[repetition][party].resize(instance.m1);
        s_prime[repetition][party].resize(instance.m1);
        t_prime[repetition][party].resize(instance.m1);
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
        for (size_t j = 0; j < instance.m1; j++) {
          vec_GF2E s_bar;
          vec_GF2E t_bar;
          // rearrange shares
          s_bar.SetLength(instance.m2 + 1);
          t_bar.SetLength(instance.m2 + 1);
          for (size_t k = 0; k < instance.m2; k++) {
            s_bar[k] = r_ejs[repetition][j] * lifted_s[j + instance.m1 * k];
            t_bar[k] = lifted_t[j + instance.m1 * k];
          }

          // sample additional random points
          std::vector<uint8_t> s_ej_bar(instance.lambda);
          random_tapes[repetition][party].squeeze_bytes(s_ej_bar.data(),
                                                        s_ej_bar.size());

          std::vector<uint8_t> t_ej_bar(instance.lambda);
          random_tapes[repetition][party].squeeze_bytes(t_ej_bar.data(),
                                                        t_ej_bar.size());
          s_bar[instance.m2] = utils::GF2E_from_bytes(s_ej_bar);
          t_bar[instance.m2] = utils::GF2E_from_bytes(t_ej_bar);

          // interpolate polynomials S_ej^i and T_ej^i
          // S_eji[repetition][party][j] =
          // utils::interpolate_with_precomputation(
          // precomputation_for_zero_to_m2, s_bar);
          // T_eji[repetition][party][j] =
          // utils::interpolate_with_precomputation(
          // precomputation_for_zero_to_m2, t_bar);
          s_prime[repetition][party][j] = s_bar;
          t_prime[repetition][party][j] = t_bar;
        }
      }
    }

    // compute sharing of P
    std::vector<vec_GF2E> P_shares(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
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
    }
    if (0 != missing_parties[repetition]) {
      for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
        // adjust first share with delta from signature
        P_shares[0][k] += proof.P_delta[k - instance.m2];
      }
    }
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // iterpolate polynomial P_e^1 from 2m+1 points
      if (party != missing_parties[repetition]) {
        P_ei[repetition][party] = utils::interpolate_with_precomputation(
            precomputation_for_zero_to_2m2, P_shares[party]);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // recompute views of polynomial checks
  /////////////////////////////////////////////////////////////////////////////
  std::vector<GF2E> c(instance.num_rounds);
  std::vector<std::vector<GF2E>> c_shares(instance.num_rounds);
  std::vector<std::vector<GF2E>> a(instance.num_rounds);
  std::vector<std::vector<std::vector<GF2E>>> a_shares(instance.num_rounds);
  std::vector<std::vector<GF2E>> b(instance.num_rounds);
  std::vector<std::vector<std::vector<GF2E>>> b_shares(instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];
    size_t missing_party = missing_parties[repetition];

    vec_GF2E lagrange_polys_evaluated_at_Re;
    lagrange_polys_evaluated_at_Re.SetLength(instance.m2 + 1);
    for (size_t k = 0; k < instance.m2 + 1; k++) {
      lagrange_polys_evaluated_at_Re[k] =
          eval(precomputation_for_zero_to_m2[k], R_es[repetition]);
    }

    c_shares[repetition].resize(instance.num_MPC_parties);
    a[repetition].resize(instance.m1);
    b[repetition].resize(instance.m1);
    a_shares[repetition].resize(instance.num_MPC_parties);
    b_shares[repetition].resize(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      a_shares[repetition][party].resize(instance.m1);
      b_shares[repetition][party].resize(instance.m1);
      if (party != missing_party) {
        for (size_t j = 0; j < instance.m1; j++) {
          // compute a_ej^i and b_ej^i
          //  a_shares[repetition][party][j] =
          //  eval(S_eji[repetition][party][j], R_es[repetition]);
          //  b_shares[repetition][party][j] =
          //  eval(T_eji[repetition][party][j], R_es[repetition]);
          a_shares[repetition][party][j] =
              lagrange_polys_evaluated_at_Re * s_prime[repetition][party][j];
          b_shares[repetition][party][j] =
              lagrange_polys_evaluated_at_Re * t_prime[repetition][party][j];
        }
        // compute c_e^i
        c_shares[repetition][party] =
            eval(P_ei[repetition][party], R_es[repetition]);
      }
    }

    // calculate missing shares
    c[repetition] = proof.P_at_R;
    c_shares[repetition][missing_party] = proof.P_at_R;
    for (size_t j = 0; j < instance.m1; j++) {
      a[repetition][j] = proof.S_j_at_R[j];
      a_shares[repetition][missing_party][j] = proof.S_j_at_R[j];
      b[repetition][j] = proof.T_j_at_R[j];
      b_shares[repetition][missing_party][j] = proof.T_j_at_R[j];
    }
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_party) {
        c_shares[repetition][missing_party] -= c_shares[repetition][party];
        for (size_t j = 0; j < instance.m1; j++) {
          a_shares[repetition][missing_party][j] -=
              a_shares[repetition][party][j];
          b_shares[repetition][missing_party][j] -=
              b_shares[repetition][party][j];
        }
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  // recompute h_1 and h_3
  /////////////////////////////////////////////////////////////////////////////
  std::vector<aes_block_t> sk_deltas;
  std::vector<std::vector<uint8_t>> t_deltas;
  for (const banquet_repetition_proof_t &proof : signature.proofs) {
    sk_deltas.push_back(proof.sk_delta);
    t_deltas.push_back(proof.t_delta);
  }
  digest_t h_1 = phase_1_commitment(instance, signature.salt, pk, message,
                                    message_len, party_seed_commitments,
                                    sk_deltas, t_deltas, rep_output_broadcasts);

  digest_t h_3 = phase_3_commitment(instance, signature.salt, h_2, c, c_shares,
                                    a, a_shares, b, b_shares);
  // do checks
  if (memcmp(h_1.data(), signature.h_1.data(), h_1.size()) != 0)
    return false;
  if (memcmp(h_2.data(), signature.h_2.data(), h_2.size()) != 0)
    return false;
  if (memcmp(h_3.data(), signature.h_3.data(), h_3.size()) != 0)
    return false;

  // check if P_e(R) = Sum_j S_e_j(R) * T_e_j(R) for all e
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    GF2E accum;
    clear(accum);
    for (size_t j = 0; j < instance.m1; j++) {
      accum += signature.proofs[repetition].S_j_at_R[j] *
               signature.proofs[repetition].T_j_at_R[j];
    }
    if (accum != signature.proofs[repetition].P_at_R)
      return false;
  }
  return true;
}

std::vector<uint8_t>
banquet_serialize_signature(const banquet_instance_t &instance,
                            const banquet_signature_t &signature) {
  std::vector<uint8_t> serialized;

  serialized.insert(serialized.end(), signature.salt.begin(),
                    signature.salt.end());
  serialized.insert(serialized.end(), signature.h_1.begin(),
                    signature.h_1.end());
  serialized.insert(serialized.end(), signature.h_2.begin(),
                    signature.h_2.end());
  serialized.insert(serialized.end(), signature.h_3.begin(),
                    signature.h_3.end());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];
    for (const seed_t &seed : proof.reveallist.first) {
      serialized.insert(serialized.end(), seed.begin(), seed.end());
    }
    serialized.insert(serialized.end(), proof.C_e.begin(), proof.C_e.end());
    serialized.insert(serialized.end(), proof.sk_delta.begin(),
                      proof.sk_delta.end());
    serialized.insert(serialized.end(), proof.t_delta.begin(),
                      proof.t_delta.end());
    serialized.insert(serialized.end(), proof.output_broadcast.begin(),
                      proof.output_broadcast.end());
    for (size_t k = 0; k < instance.m2 + 1; k++) {
      std::vector<uint8_t> buffer = GF2E_to_bytes(instance, proof.P_delta[k]);
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
    {
      std::vector<uint8_t> buffer = GF2E_to_bytes(instance, proof.P_at_R);
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer = GF2E_to_bytes(instance, proof.S_j_at_R[j]);
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer = GF2E_to_bytes(instance, proof.T_j_at_R[j]);
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
  }
  return serialized;
}
banquet_signature_t
banquet_deserialize_signature(const banquet_instance_t &instance,
                              const std::vector<uint8_t> &serialized) {

  size_t current_offset = 0;
  banquet_salt_t salt;
  memcpy(salt.data(), serialized.data() + current_offset, salt.size());
  current_offset += salt.size();
  digest_t h_1, h_2, h_3;
  memcpy(h_1.data(), serialized.data() + current_offset, h_1.size());
  current_offset += h_1.size();
  memcpy(h_2.data(), serialized.data() + current_offset, h_2.size());
  current_offset += h_2.size();
  memcpy(h_3.data(), serialized.data() + current_offset, h_3.size());
  current_offset += h_3.size();
  std::vector<banquet_repetition_proof_t> proofs;

  std::vector<uint8_t> missing_parties = phase_3_expand(instance, h_3);
  size_t reveallist_size = ceil_log2(instance.num_MPC_parties);
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    reveal_list_t reveallist;
    reveallist.first.reserve(reveallist_size);
    reveallist.second = missing_parties[repetition];
    for (size_t i = 0; i < reveallist_size; i++) {
      seed_t seed;
      memcpy(seed.data(), serialized.data() + current_offset, seed.size());
      current_offset += seed.size();
      reveallist.first.push_back(seed);
    }
    digest_t C_e;
    memcpy(C_e.data(), serialized.data() + current_offset, C_e.size());
    current_offset += C_e.size();

    aes_block_t sk_delta;
    memcpy(sk_delta.data(), serialized.data() + current_offset,
           sk_delta.size());
    current_offset += sk_delta.size();

    std::vector<uint8_t> t_delta(NUM_SBOXES_AES_128);
    memcpy(t_delta.data(), serialized.data() + current_offset, t_delta.size());
    current_offset += t_delta.size();

    aes_block_t output_broadcast;
    memcpy(output_broadcast.data(), serialized.data() + current_offset,
           output_broadcast.size());
    current_offset += output_broadcast.size();

    std::vector<GF2E> P_delta;
    P_delta.reserve(instance.m2 + 1);
    for (size_t k = 0; k < instance.m2 + 1; k++) {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      P_delta.push_back(utils::GF2E_from_bytes(buffer));
    }
    GF2E P_at_R;
    {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      P_at_R = utils::GF2E_from_bytes(buffer);
    }
    std::vector<GF2E> S_j_at_R;
    S_j_at_R.reserve(instance.m1);
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      S_j_at_R.push_back(utils::GF2E_from_bytes(buffer));
    }
    std::vector<GF2E> T_j_at_R;
    T_j_at_R.reserve(instance.m1);
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      T_j_at_R.push_back(utils::GF2E_from_bytes(buffer));
    }
    proofs.emplace_back(banquet_repetition_proof_t{
        reveallist, C_e, sk_delta, t_delta, output_broadcast, P_delta, P_at_R,
        S_j_at_R, T_j_at_R});
  }
  assert(current_offset == serialized.size());
  banquet_signature_t signature{salt, h_1, h_2, h_3, proofs};
  return signature;
}