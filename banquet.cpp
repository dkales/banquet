#include "banquet.h"

#include "aes.h"
#include "field.h"
#include "tape.h"
#include "tree.h"
#include <algorithm>
#include <cassert>
#include <cstring>

extern "C" {
#include "kdf_shake.h"
#include "randomness.h"
}

namespace {
inline void hash_update_GF2E(hash_context *ctx,
                             const banquet_instance_t &instance,
                             const field::GF2E &element) {
  // 8 bytes is enough for supported field sizes
  std::array<uint8_t, 8> buffer;
  element.to_bytes(buffer.data());
  hash_update(ctx, buffer.data(), instance.lambda);
}

std::pair<banquet_salt_t, std::vector<std::vector<uint8_t>>>
generate_salt_and_seeds(const banquet_instance_t &instance,
                        const banquet_keypair_t &keypair,
                        const uint8_t *message, size_t message_len) {
  // salt, seed_1, ..., seed_r = H(instance||sk||pk||m)
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update_uint16_le(&ctx, (uint16_t)instance.params);
  hash_update(&ctx, keypair.first.data(), keypair.first.size());
  hash_update(&ctx, keypair.second.data(), keypair.second.size());
  hash_update(&ctx, message, message_len);
  hash_final(&ctx);

  banquet_salt_t salt;
  hash_squeeze(&ctx, salt.data(), salt.size());
  std::vector<std::vector<uint8_t>> seeds;
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    std::vector<uint8_t> s(instance.seed_size);
    hash_squeeze(&ctx, s.data(), s.size());
    seeds.push_back(s);
  }
  return std::make_pair(salt, seeds);
}

std::vector<uint8_t> commit_to_party_seed(const banquet_instance_t &instance,
                                          const std::vector<uint8_t> &seed,
                                          const banquet_salt_t &salt,
                                          size_t rep_idx, size_t party_idx) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_idx);
  hash_update_uint16_le(&ctx, (uint16_t)party_idx);
  hash_update(&ctx, seed.data(), seed.size());
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<uint8_t> phase_1_commitment(
    const banquet_instance_t &instance, const banquet_salt_t &salt,
    const std::vector<uint8_t> &pk, const uint8_t *message, size_t message_len,
    const std::vector<std::vector<std::vector<uint8_t>>> &commitments,
    const std::vector<std::vector<uint8_t>> &key_deltas,
    const std::vector<std::vector<uint8_t>> &t_deltas,
    RepByteContainer &output_broadcasts) {

  hash_context ctx;
  hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_1);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, pk.data(), pk.size());
  hash_update(&ctx, message, message_len);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      hash_update(&ctx, commitments[repetition][party].data(),
                  commitments[repetition][party].size());
      auto output_broadcast = output_broadcasts.get(repetition, party);
      hash_update(&ctx, output_broadcast.data(), output_broadcast.size());
    }
    hash_update(&ctx, key_deltas[repetition].data(),
                key_deltas[repetition].size());
    hash_update(&ctx, t_deltas[repetition].data(), t_deltas[repetition].size());
  }
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<std::vector<field::GF2E>>
phase_1_expand(const banquet_instance_t &instance,
               const std::vector<uint8_t> &h_1) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_1.data(), h_1.size());
  hash_final(&ctx);

  std::vector<std::vector<field::GF2E>> r_ejs;
  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);
  r_ejs.reserve(instance.num_rounds);
  for (size_t e = 0; e < instance.num_rounds; e++) {
    std::vector<field::GF2E> r_js;
    r_js.resize(instance.m1);
    for (size_t j = 0; j < instance.m1; j++) {
      hash_squeeze(&ctx, lambda_sized_buffer.data(),
                   lambda_sized_buffer.size());
      r_js[j].from_bytes(lambda_sized_buffer.data());
    }
    r_ejs.push_back(r_js);
  }
  return r_ejs;
}

std::vector<uint8_t>
phase_2_commitment(const banquet_instance_t &instance,
                   const banquet_salt_t &salt, const std::vector<uint8_t> &h_1,
                   const std::vector<std::vector<field::GF2E>> &P_deltas) {

  hash_context ctx;
  hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_2);
  hash_update(&ctx, salt.data(), salt.size());
  hash_update(&ctx, h_1.data(), h_1.size());

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    for (size_t k = 0; k < instance.m2; k++) {
      hash_update_GF2E(&ctx, instance, P_deltas[repetition][k]);
    }
  }
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<field::GF2E>
phase_2_expand(const banquet_instance_t &instance,
               const std::vector<uint8_t> &h_2,
               const std::vector<field::GF2E> &forbidden_values) {
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
  hash_update(&ctx, h_2.data(), h_2.size());
  hash_final(&ctx);

  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);
  std::vector<field::GF2E> R_es;
  for (size_t e = 0; e < instance.num_rounds; e++) {
    while (true) {
      hash_squeeze(&ctx, lambda_sized_buffer.data(),
                   lambda_sized_buffer.size());
      //  check that R is not in {0,...m2-1}
      field::GF2E candidate_R;
      candidate_R.from_bytes(lambda_sized_buffer.data());
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

std::vector<uint8_t>
phase_3_commitment(const banquet_instance_t &instance,
                   const banquet_salt_t &salt, const std::vector<uint8_t> &h_2,
                   const std::vector<field::GF2E> &c,
                   const std::vector<std::vector<field::GF2E>> &c_shares,
                   const std::vector<std::vector<field::GF2E>> &a,
                   const RepContainer<field::GF2E> &a_shares,
                   const std::vector<std::vector<field::GF2E>> &b,
                   const RepContainer<field::GF2E> &b_shares) {

  hash_context ctx;
  hash_init_prefix(&ctx, instance.digest_size, HASH_PREFIX_3);
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
        hash_update_GF2E(&ctx, instance, a_shares.get(repetition, party)[j]);
        hash_update_GF2E(&ctx, instance, b_shares.get(repetition, party)[j]);
      }
    }
  }
  hash_final(&ctx);

  std::vector<uint8_t> commitment(instance.digest_size);
  hash_squeeze(&ctx, commitment.data(), commitment.size());
  return commitment;
}

std::vector<uint8_t> phase_3_expand(const banquet_instance_t &instance,
                                    const std::vector<uint8_t> &h_3) {
  assert(instance.num_MPC_parties <= 256);
  // TODO assert(is_power_of_2(instance.num_MPC_parties));
  hash_context ctx;
  hash_init(&ctx, instance.digest_size);
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
  std::vector<uint8_t> key(instance.aes_params.key_size),
      pt(instance.aes_params.block_size * instance.aes_params.num_blocks),
      ct(instance.aes_params.block_size * instance.aes_params.num_blocks);

  while (true) {
    rand_bytes(key.data(), key.size());
    rand_bytes(pt.data(), pt.size());
    if (instance.aes_params.key_size == 16) {
      if (AES128::aes_128(key, pt, ct)) {
        break;
      }
    } else if (instance.aes_params.key_size == 24) {
      if (AES192::aes_192(key, pt, ct)) {
        break;
      }
    } else if (instance.aes_params.key_size == 32) {
      if (AES256::aes_256(key, pt, ct)) {
        break;
      }
    } else
      throw std::runtime_error("invalid parameters");
  }
  banquet_keypair_t keypair;
  keypair.first = key;
  keypair.second = pt;
  keypair.second.insert(keypair.second.end(), ct.begin(), ct.end());
  return keypair;
}

banquet_signature_t banquet_sign(const banquet_instance_t &instance,
                                 const banquet_keypair_t &keypair,
                                 const uint8_t *message, size_t message_len) {
  // init modulus of extension field F_{2^{8\lambda}}
  field::GF2E::init_extension_field(instance);

  // grab aes key, pt and ct
  std::vector<uint8_t> key = keypair.first;
  std::vector<uint8_t> pt_ct = keypair.second;
  const size_t total_pt_ct_size =
      instance.aes_params.block_size * instance.aes_params.num_blocks;
  std::vector<uint8_t> pt(total_pt_ct_size), ct(total_pt_ct_size),
      ct2(total_pt_ct_size);
  memcpy(pt.data(), keypair.second.data(), pt.size());
  memcpy(ct.data(), keypair.second.data() + pt.size(), ct.size());

  // get sbox inputs and outputs for aes evaluation
  std::pair<std::vector<uint8_t>, std::vector<uint8_t>> sbox_pairs;

  if (instance.aes_params.key_size == 16)
    sbox_pairs = AES128::aes_128_with_sbox_output(key, pt, ct2);
  else if (instance.aes_params.key_size == 24)
    sbox_pairs = AES192::aes_192_with_sbox_output(key, pt, ct2);
  else if (instance.aes_params.key_size == 32)
    sbox_pairs = AES256::aes_256_with_sbox_output(key, pt, ct2);
  else
    throw std::runtime_error("invalid parameters");
  // sanity check, incoming keypair is valid
  assert(ct == ct2);

  // generate salt and master seeds for each repetition
  auto [salt, master_seeds] =
      generate_salt_and_seeds(instance, keypair, message, message_len);

  // buffer for squeezing field elements into
  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);

  // do parallel repetitions
  // create seed trees and random tapes
  std::vector<SeedTree> seed_trees;
  std::vector<std::vector<RandomTape>> random_tapes;
  std::vector<std::vector<std::vector<uint8_t>>> party_seed_commitments;

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    // generate seed tree for the N parties
    seed_trees.emplace_back(master_seeds[repetition], instance.num_MPC_parties,
                            salt, repetition);

    // commit to each party's seed;
    std::vector<std::vector<uint8_t>> current_party_seed_commitments;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      current_party_seed_commitments.push_back(commit_to_party_seed(
          instance, seed_trees[repetition].get_leaf(party).value(), salt,
          repetition, party));
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
  RepByteContainer rep_shared_keys(instance.num_rounds,
                                   instance.num_MPC_parties,
                                   instance.aes_params.key_size);
  RepByteContainer rep_output_broadcasts(
      instance.num_rounds, instance.num_MPC_parties,
      instance.aes_params.block_size * instance.aes_params.num_blocks);
  RepByteContainer rep_shared_s(instance.num_rounds, instance.num_MPC_parties,
                                instance.aes_params.num_sboxes);
  RepByteContainer rep_shared_t(instance.num_rounds, instance.num_MPC_parties,
                                instance.aes_params.num_sboxes);
  std::vector<std::vector<uint8_t>> rep_key_deltas;
  std::vector<std::vector<uint8_t>> rep_t_deltas;

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {

    // generate sharing of secret key
    std::vector<uint8_t> key_delta = key;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      auto shared_key = rep_shared_keys.get(repetition, party);
      random_tapes[repetition][party].squeeze_bytes(shared_key.data(),
                                                    shared_key.size());
      std::transform(std::begin(shared_key), std::end(shared_key),
                     std::begin(key_delta), std::begin(key_delta),
                     std::bit_xor<uint8_t>());
    }

    // fix first share
    auto first_share_key = rep_shared_keys.get(repetition, 0);
    std::transform(std::begin(key_delta), std::end(key_delta),
                   std::begin(first_share_key), std::begin(first_share_key),
                   std::bit_xor<uint8_t>());

    rep_key_deltas.push_back(key_delta);
    // generate sharing of t values
    std::vector<uint8_t> t_deltas = sbox_pairs.second;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      auto shared_t = rep_shared_t.get(repetition, party);
      random_tapes[repetition][party].squeeze_bytes(shared_t.data(),
                                                    shared_t.size());
      std::transform(std::begin(shared_t), std::end(shared_t),
                     std::begin(t_deltas), std::begin(t_deltas),
                     std::bit_xor<uint8_t>());
    }
    // fix first share
    auto first_share_t = rep_shared_t.get(repetition, 0);
    std::transform(std::begin(t_deltas), std::end(t_deltas),
                   std::begin(first_share_t), std::begin(first_share_t),
                   std::bit_xor<uint8_t>());

    // get shares of sbox inputs by executing MPC AES
    auto ct_shares = rep_output_broadcasts.get_repetition(repetition);
    auto shared_s = rep_shared_s.get_repetition(repetition);

    if (instance.aes_params.key_size == 16)
      AES128::aes_128_s_shares(rep_shared_keys.get_repetition(repetition),
                               rep_shared_t.get_repetition(repetition), pt,
                               ct_shares, shared_s);
    else if (instance.aes_params.key_size == 24)
      AES192::aes_192_s_shares(rep_shared_keys.get_repetition(repetition),
                               rep_shared_t.get_repetition(repetition), pt,
                               ct_shares, shared_s);
    else if (instance.aes_params.key_size == 32)
      AES256::aes_256_s_shares(rep_shared_keys.get_repetition(repetition),
                               rep_shared_t.get_repetition(repetition), pt,
                               ct_shares, shared_s);
    else
      throw std::runtime_error("invalid parameters");

    // sanity check, mpc execution = plain one
    std::vector<uint8_t> ct_check(instance.aes_params.block_size *
                                  instance.aes_params.num_blocks);
    memset(ct_check.data(), 0, ct_check.size());
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      std::transform(std::begin(ct_shares[party]), std::end(ct_shares[party]),
                     std::begin(ct_check), std::begin(ct_check),
                     std::bit_xor<uint8_t>());
    }

    assert(ct == ct_check);
    rep_t_deltas.push_back(t_deltas);
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 2: challenge the multiplications
  /////////////////////////////////////////////////////////////////////////////

  // commit to salt, (all commitments of parties seeds, key_delta, t_delta)
  // for all repetitions
  std::vector<uint8_t> h_1 =
      phase_1_commitment(instance, salt, keypair.second, message, message_len,
                         party_seed_commitments, rep_key_deltas, rep_t_deltas,
                         rep_output_broadcasts);

  // expand challenge hash to M * m1 values
  std::vector<std::vector<field::GF2E>> r_ejs = phase_1_expand(instance, h_1);

  /////////////////////////////////////////////////////////////////////////////
  // phase 3: commit to the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  // a vector of the first m2+1 field elements for interpolation
  std::vector<field::GF2E> x_values_for_interpolation_zero_to_m2 =
      field::get_first_n_field_elements(instance.m2 + 1);
  std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_m2 =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_m2);
  std::vector<field::GF2E> x_values_for_interpolation_zero_to_2m2 =
      field::get_first_n_field_elements(2 * instance.m2 + 1);
  std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2m2 =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_2m2);

  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> s_prime(
      instance.num_rounds);
  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> t_prime(
      instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> S_eji(instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> T_eji(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> P_e(instance.num_rounds);
  std::vector<std::vector<std::vector<field::GF2E>>> P_e_shares(
      instance.num_rounds);
  // std::vector<std::vector<GF2EX>> P_ei(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> P_deltas(instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    s_prime[repetition].resize(instance.num_MPC_parties);
    t_prime[repetition].resize(instance.num_MPC_parties);
    // S_eji[repetition].resize(instance.num_MPC_parties);
    // T_eji[repetition].resize(instance.num_MPC_parties);
    // P_ei[repetition].resize(instance.num_MPC_parties);
    P_deltas[repetition].resize(instance.m2 + 1);

    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      s_prime[repetition][party].resize(instance.m1);
      t_prime[repetition][party].resize(instance.m1);
      // S_eji[repetition][party].resize(instance.m1);
      // T_eji[repetition][party].resize(instance.m1);
      // lift shares from F_{2^8} to F_{2^{8\lambda}}
      std::vector<std::reference_wrapper<const field::GF2E>> lifted_s;
      std::vector<std::reference_wrapper<const field::GF2E>> lifted_t;
      lifted_s.reserve(instance.aes_params.num_sboxes);
      lifted_t.reserve(instance.aes_params.num_sboxes);
      auto shared_s = rep_shared_s.get(repetition, party);
      auto shared_t = rep_shared_t.get(repetition, party);
      for (size_t idx = 0; idx < instance.aes_params.num_sboxes; idx++) {
        lifted_s.push_back(field::lift_uint8_t(shared_s[idx]));
        lifted_t.push_back(field::lift_uint8_t(shared_t[idx]));
      }

      // rearrange shares
      std::vector<field::GF2E> s_bar;
      std::vector<field::GF2E> t_bar;
      s_bar.resize(instance.m2 + 1);
      t_bar.resize(instance.m2 + 1);

      for (size_t j = 0; j < instance.m1; j++) {
        for (size_t k = 0; k < instance.m2; k++) {
          s_bar[k] = r_ejs[repetition][j] * lifted_s[j + instance.m1 * k];
          t_bar[k] = lifted_t[j + instance.m1 * k];
        }

        // sample additional random points
        random_tapes[repetition][party].squeeze_bytes(
            lambda_sized_buffer.data(), lambda_sized_buffer.size());
        s_bar[instance.m2].from_bytes(lambda_sized_buffer.data());

        random_tapes[repetition][party].squeeze_bytes(
            lambda_sized_buffer.data(), lambda_sized_buffer.size());
        t_bar[instance.m2].from_bytes(lambda_sized_buffer.data());

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
    std::vector<field::GF2E> P(2 * instance.m2 + 1);
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<field::GF2E> s_prime_sum(instance.m2 + 1);
      std::vector<field::GF2E> t_prime_sum(instance.m2 + 1);
      std::vector<field::GF2E> S_sum;
      std::vector<field::GF2E> T_sum;

      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        // S_sum += S_eji[repetition][party][j];
        // T_sum += T_eji[repetition][party][j];
        s_prime_sum += s_prime[repetition][party][j];
        t_prime_sum += t_prime[repetition][party][j];
      }
      S_sum = field::interpolate_with_precomputation(
          precomputation_for_zero_to_m2, s_prime_sum);
      T_sum = field::interpolate_with_precomputation(
          precomputation_for_zero_to_m2, t_prime_sum);

      P += S_sum * T_sum;
    }
    P_e[repetition] = P;

    // compute sharing of P
    std::vector<std::vector<field::GF2E>> &P_shares = P_e_shares[repetition];
    P_shares.resize(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      // first m2 points: first party = sum of r_e,j, other parties = 0
      P_shares[party].resize(2 * instance.m2 + 1);
      if (party == 0) {
        field::GF2E sum_r;
        for (size_t j = 0; j < instance.m1; j++) {
          sum_r += r_ejs[repetition][j];
        }
        for (size_t k = 0; k < instance.m2; k++) {
          P_shares[party][k] = sum_r;
        }
      } else {
        for (size_t k = 0; k < instance.m2; k++) {
          P_shares[party][k] = field::GF2E(0);
        }
      }

      // second m2+1 points: sample from random tape
      for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
        random_tapes[repetition][party].squeeze_bytes(
            lambda_sized_buffer.data(), lambda_sized_buffer.size());
        P_shares[party][k].from_bytes(lambda_sized_buffer.data());
      }
    }
    for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
      // calculate offset
      field::GF2E k_element = x_values_for_interpolation_zero_to_2m2[k];
      field::GF2E P_at_k_delta = field::eval(P, k_element);
      for (size_t party = 0; party < instance.num_MPC_parties; party++) {
        P_at_k_delta -= P_shares[party][k];
      }
      P_deltas[repetition][k - instance.m2] = P_at_k_delta;
      // adjust first share
      P_shares[0][k] += P_at_k_delta;
    }
    // for (size_t party = 0; party < instance.num_MPC_parties; party++) {
    //// iterpolate polynomial P_e^1 from 2m+1 points
    // P_ei[repetition][party] = utils::interpolate_with_precomputation(
    // precomputation_for_zero_to_2m2, P_shares[party]);
    //}
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 4: challenge the checking polynomials
  /////////////////////////////////////////////////////////////////////////////

  std::vector<uint8_t> h_2 = phase_2_commitment(instance, salt, h_1, P_deltas);

  // expand challenge hash to M values

  std::vector<field::GF2E> forbidden_challenge_values =
      field::get_first_n_field_elements(instance.m2);
  std::vector<field::GF2E> R_es =
      phase_2_expand(instance, h_2, forbidden_challenge_values);

  /////////////////////////////////////////////////////////////////////////////
  // phase 5: commit to the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  std::vector<field::GF2E> c(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> c_shares(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> a(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> b(instance.num_rounds);
  RepContainer<field::GF2E> a_shares(instance.num_rounds,
                                     instance.num_MPC_parties, instance.m1);
  RepContainer<field::GF2E> b_shares(instance.num_rounds,
                                     instance.num_MPC_parties, instance.m1);

  std::vector<field::GF2E> lagrange_polys_evaluated_at_Re_m2(instance.m2 + 1);
  std::vector<field::GF2E> lagrange_polys_evaluated_at_Re_2m2(2 * instance.m2 +
                                                              1);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    c_shares[repetition].resize(instance.num_MPC_parties);
    a[repetition].resize(instance.m1);
    b[repetition].resize(instance.m1);
    for (size_t k = 0; k < instance.m2 + 1; k++) {
      lagrange_polys_evaluated_at_Re_m2[k] =
          field::eval(precomputation_for_zero_to_m2[k], R_es[repetition]);
    }
    for (size_t k = 0; k < 2 * instance.m2 + 1; k++) {
      lagrange_polys_evaluated_at_Re_2m2[k] =
          field::eval(precomputation_for_zero_to_2m2[k], R_es[repetition]);
    }

    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      auto a_shares_party = a_shares.get(repetition, party);
      auto b_shares_party = b_shares.get(repetition, party);
      for (size_t j = 0; j < instance.m1; j++) {
        // compute a_ej^i and b_ej^i
        // a_shares[repetition][party][j] =
        // eval(S_eji[repetition][party][j], R_es[repetition]);
        // b_shares[repetition][party][j] =
        // eval(T_eji[repetition][party][j], R_es[repetition]);
        a_shares_party[j] = dot_product(lagrange_polys_evaluated_at_Re_m2,
                                        s_prime[repetition][party][j]);
        b_shares_party[j] = dot_product(lagrange_polys_evaluated_at_Re_m2,
                                        t_prime[repetition][party][j]);
      }
      // compute c_e^i
      c_shares[repetition][party] = dot_product(
          lagrange_polys_evaluated_at_Re_2m2, P_e_shares[repetition][party]);
    }
    // open c_e and a,b values
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      c[repetition] += c_shares[repetition][party];
      auto a_shares_party = a_shares.get(repetition, party);
      auto b_shares_party = b_shares.get(repetition, party);
      for (size_t j = 0; j < instance.m1; j++) {
        a[repetition][j] += a_shares_party[j];
        b[repetition][j] += b_shares_party[j];
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // phase 6: challenge the views of the checking protocol
  /////////////////////////////////////////////////////////////////////////////

  std::vector<uint8_t> h_3 = phase_3_commitment(
      instance, salt, h_2, c, c_shares, a, a_shares, b, b_shares);

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
        std::vector<uint8_t>(
            rep_output_broadcasts.get(repetition, missing_party).begin(),
            rep_output_broadcasts.get(repetition, missing_party).end()),
        P_deltas[repetition],
        c[repetition],
        a[repetition],
        b[repetition],
    };
    // sanity check c = sum_j a*b
    field::GF2E accum;
    for (size_t j = 0; j < instance.m1; j++) {
      accum += a[repetition][j] * b[repetition][j];
    }
    if (accum != c[repetition])
      throw std::runtime_error("something wrong here");
    proofs.push_back(proof);
  }

  banquet_signature_t signature{salt, h_1, h_2, h_3, proofs};

  return signature;
}

bool banquet_verify(const banquet_instance_t &instance,
                    const std::vector<uint8_t> &pk,
                    const banquet_signature_t &signature,
                    const uint8_t *message, size_t message_len) {

  // init modulus of extension field F_{2^{8\lambda}}
  field::GF2E::init_extension_field(instance);

  std::vector<uint8_t> pt(instance.aes_params.block_size *
                          instance.aes_params.num_blocks),
      ct(instance.aes_params.block_size * instance.aes_params.num_blocks);
  memcpy(pt.data(), pk.data(), pt.size());
  memcpy(ct.data(), pk.data() + pt.size(), ct.size());

  // buffer for squeezing field elements into
  std::vector<uint8_t> lambda_sized_buffer(instance.lambda);

  // do parallel repetitions
  // create seed trees and random tapes
  std::vector<SeedTree> seed_trees;
  std::vector<std::vector<RandomTape>> random_tapes;
  std::vector<std::vector<std::vector<uint8_t>>> party_seed_commitments;

  // recompute h_2
  std::vector<std::vector<field::GF2E>> P_deltas;
  for (const banquet_repetition_proof_t &proof : signature.proofs) {
    P_deltas.push_back(proof.P_delta);
  }
  std::vector<uint8_t> h_2 =
      phase_2_commitment(instance, signature.salt, signature.h_1, P_deltas);

  // compute challenges based on hashes
  // h1 expansion
  std::vector<std::vector<field::GF2E>> r_ejs =
      phase_1_expand(instance, signature.h_1);
  // h2 expansion
  std::vector<field::GF2E> forbidden_challenge_values =
      field::get_first_n_field_elements(instance.m2);
  std::vector<field::GF2E> R_es =
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
    std::vector<std::vector<uint8_t>> current_party_seed_commitments;
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        current_party_seed_commitments.push_back(commit_to_party_seed(
            instance, seed_trees[repetition].get_leaf(party).value(),
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
        party_tapes.emplace_back(std::vector<uint8_t>(instance.seed_size),
                                 signature.salt, repetition, party);
      }
    }
    random_tapes.push_back(party_tapes);
  }
  /////////////////////////////////////////////////////////////////////////////
  // recompute commitments to executions of AES
  /////////////////////////////////////////////////////////////////////////////
  RepByteContainer rep_shared_keys(instance.num_rounds,
                                   instance.num_MPC_parties,
                                   instance.aes_params.key_size);
  RepByteContainer rep_shared_s(instance.num_rounds, instance.num_MPC_parties,
                                instance.aes_params.num_sboxes);
  RepByteContainer rep_shared_t(instance.num_rounds, instance.num_MPC_parties,
                                instance.aes_params.num_sboxes);
  RepByteContainer rep_output_broadcasts(
      instance.num_rounds, instance.num_MPC_parties,
      instance.aes_params.block_size * instance.aes_params.num_blocks);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];

    // generate sharing of secret key
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      auto shared_key = rep_shared_keys.get(repetition, party);
      random_tapes[repetition][party].squeeze_bytes(shared_key.data(),
                                                    shared_key.size());
    }

    // fix first share
    auto first_key_share = rep_shared_keys.get(repetition, 0);
    std::transform(std::begin(proof.sk_delta), std::end(proof.sk_delta),
                   std::begin(first_key_share), std::begin(first_key_share),
                   std::bit_xor<uint8_t>());

    // generate sharing of t values
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      auto shared_t = rep_shared_t.get(repetition, party);
      random_tapes[repetition][party].squeeze_bytes(shared_t.data(),
                                                    shared_t.size());
    }
    // fix first share
    auto first_shared_t = rep_shared_t.get(repetition, 0);
    std::transform(std::begin(proof.t_delta), std::end(proof.t_delta),
                   std::begin(first_shared_t), std::begin(first_shared_t),
                   std::bit_xor<uint8_t>());

    // get shares of sbox inputs by executing MPC AES
    auto ct_shares = rep_output_broadcasts.get_repetition(repetition);
    auto shared_s = rep_shared_s.get_repetition(repetition);

    if (instance.aes_params.key_size == 16)
      AES128::aes_128_s_shares(rep_shared_keys.get_repetition(repetition),
                               rep_shared_t.get_repetition(repetition), pt,
                               ct_shares, shared_s);
    else if (instance.aes_params.key_size == 24)
      AES192::aes_192_s_shares(rep_shared_keys.get_repetition(repetition),
                               rep_shared_t.get_repetition(repetition), pt,
                               ct_shares, shared_s);
    else if (instance.aes_params.key_size == 32)
      AES256::aes_256_s_shares(rep_shared_keys.get_repetition(repetition),
                               rep_shared_t.get_repetition(repetition), pt,
                               ct_shares, shared_s);
    else
      throw std::runtime_error("invalid parameters");

    // get missing output broadcast from proof
    std::copy(proof.output_broadcast.begin(), proof.output_broadcast.end(),
              ct_shares[missing_parties[repetition]].begin());
    // check MPC execution is correct
    std::vector<uint8_t> ct_check(instance.aes_params.block_size *
                                  instance.aes_params.num_blocks);
    memset(ct_check.data(), 0, ct_check.size());
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      std::transform(std::begin(ct_shares[party]), std::end(ct_shares[party]),
                     std::begin(ct_check), std::begin(ct_check),
                     std::bit_xor<uint8_t>());
    }
    if (memcmp(ct_check.data(), ct.data(), ct.size()) != 0) {
      return false;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // recompute shares of polynomials
  /////////////////////////////////////////////////////////////////////////////
  // a vector of the first m2+1 field elements for interpolation
  std::vector<field::GF2E> x_values_for_interpolation_zero_to_m2 =
      field::get_first_n_field_elements(instance.m2 + 1);
  std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_m2 =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_m2);
  std::vector<field::GF2E> x_values_for_interpolation_zero_to_2m2 =
      field::get_first_n_field_elements(2 * instance.m2 + 1);
  std::vector<std::vector<field::GF2E>> precomputation_for_zero_to_2m2 =
      field::precompute_lagrange_polynomials(
          x_values_for_interpolation_zero_to_2m2);

  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> s_prime(
      instance.num_rounds);
  std::vector<std::vector<std::vector<std::vector<field::GF2E>>>> t_prime(
      instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> S_eji(instance.num_rounds);
  // std::vector<std::vector<std::vector<GF2EX>>> T_eji(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> P_e(instance.num_rounds);
  // std::vector<std::vector<GF2EX>> P_ei(instance.num_rounds);
  std::vector<std::vector<std::vector<field::GF2E>>> P_e_shares(
      instance.num_rounds);

  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];
    // S_eji[repetition].resize(instance.num_MPC_parties);
    // T_eji[repetition].resize(instance.num_MPC_parties);
    s_prime[repetition].resize(instance.num_MPC_parties);
    t_prime[repetition].resize(instance.num_MPC_parties);
    // P_ei[repetition].resize(instance.num_MPC_parties);
    P_deltas[repetition].resize(instance.m2 + 1);

    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        // S_eji[repetition][party].resize(instance.m1);
        // T_eji[repetition][party].resize(instance.m1);
        s_prime[repetition][party].resize(instance.m1);
        t_prime[repetition][party].resize(instance.m1);
        // lift shares from F_{2^8} to F_{2^{8\lambda}}
        std::vector<std::reference_wrapper<const field::GF2E>> lifted_s;
        std::vector<std::reference_wrapper<const field::GF2E>> lifted_t;
        lifted_s.reserve(instance.aes_params.num_sboxes);
        lifted_t.reserve(instance.aes_params.num_sboxes);
        auto shared_s = rep_shared_s.get(repetition, party);
        auto shared_t = rep_shared_t.get(repetition, party);
        for (size_t idx = 0; idx < instance.aes_params.num_sboxes; idx++) {
          lifted_s.push_back(field::lift_uint8_t(shared_s[idx]));
          lifted_t.push_back(field::lift_uint8_t(shared_t[idx]));
        }

        // rearrange shares
        std::vector<field::GF2E> s_bar(instance.m2 + 1);
        std::vector<field::GF2E> t_bar(instance.m2 + 1);

        for (size_t j = 0; j < instance.m1; j++) {
          for (size_t k = 0; k < instance.m2; k++) {
            s_bar[k] = r_ejs[repetition][j] * lifted_s[j + instance.m1 * k];
            t_bar[k] = lifted_t[j + instance.m1 * k];
          }

          // sample additional random points
          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          s_bar[instance.m2].from_bytes(lambda_sized_buffer.data());

          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          t_bar[instance.m2].from_bytes(lambda_sized_buffer.data());

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
    std::vector<std::vector<field::GF2E>> &P_shares = P_e_shares[repetition];
    P_shares.resize(instance.num_MPC_parties);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_parties[repetition]) {
        // first m2 points: first party = sum of r_e,j, other parties = 0
        P_shares[party].resize(2 * instance.m2 + 1);
        if (party == 0) {
          field::GF2E sum_r;
          for (size_t j = 0; j < instance.m1; j++) {
            sum_r += r_ejs[repetition][j];
          }
          for (size_t k = 0; k < instance.m2; k++) {
            P_shares[party][k] = sum_r;
          }
        } else {
          for (size_t k = 0; k < instance.m2; k++) {
            P_shares[party][k] = field::GF2E(0);
          }
        }

        // second m2+1 points: sample from random tape
        for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
          random_tapes[repetition][party].squeeze_bytes(
              lambda_sized_buffer.data(), lambda_sized_buffer.size());
          P_shares[party][k].from_bytes(lambda_sized_buffer.data());
        }
      }
    }
    if (0 != missing_parties[repetition]) {
      for (size_t k = instance.m2; k <= 2 * instance.m2; k++) {
        // adjust first share with delta from signature
        P_shares[0][k] += proof.P_delta[k - instance.m2];
      }
    }
    // for (size_t party = 0; party < instance.num_MPC_parties; party++) {
    //// iterpolate polynomial P_e^1 from 2m+1 points
    // if (party != missing_parties[repetition]) {
    // P_ei[repetition][party] = utils::interpolate_with_precomputation(
    // precomputation_for_zero_to_2m2, P_shares[party]);
    //}
    //}
  }

  /////////////////////////////////////////////////////////////////////////////
  // recompute views of polynomial checks
  /////////////////////////////////////////////////////////////////////////////
  std::vector<field::GF2E> c(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> c_shares(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> a(instance.num_rounds);
  std::vector<std::vector<field::GF2E>> b(instance.num_rounds);
  RepContainer<field::GF2E> a_shares(instance.num_rounds,
                                     instance.num_MPC_parties, instance.m1);
  RepContainer<field::GF2E> b_shares(instance.num_rounds,
                                     instance.num_MPC_parties, instance.m1);

  std::vector<field::GF2E> lagrange_polys_evaluated_at_Re_m2(instance.m2 + 1);
  std::vector<field::GF2E> lagrange_polys_evaluated_at_Re_2m2(2 * instance.m2 +
                                                              1);
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    const banquet_repetition_proof_t &proof = signature.proofs[repetition];
    size_t missing_party = missing_parties[repetition];

    for (size_t k = 0; k < instance.m2 + 1; k++) {
      lagrange_polys_evaluated_at_Re_m2[k] =
          field::eval(precomputation_for_zero_to_m2[k], R_es[repetition]);
    }
    for (size_t k = 0; k < 2 * instance.m2 + 1; k++) {
      lagrange_polys_evaluated_at_Re_2m2[k] =
          field::eval(precomputation_for_zero_to_2m2[k], R_es[repetition]);
    }

    c_shares[repetition].resize(instance.num_MPC_parties);
    a[repetition].resize(instance.m1);
    b[repetition].resize(instance.m1);
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_party) {
        auto a_shares_party = a_shares.get(repetition, party);
        auto b_shares_party = b_shares.get(repetition, party);
        for (size_t j = 0; j < instance.m1; j++) {
          // compute a_ej^i and b_ej^i
          //  a_shares[repetition][party][j] =
          //  eval(S_eji[repetition][party][j], R_es[repetition]);
          //  b_shares[repetition][party][j] =
          //  eval(T_eji[repetition][party][j], R_es[repetition]);
          a_shares_party[j] = dot_product(lagrange_polys_evaluated_at_Re_m2,
                                          s_prime[repetition][party][j]);
          b_shares_party[j] = dot_product(lagrange_polys_evaluated_at_Re_m2,
                                          t_prime[repetition][party][j]);
        }
        // compute c_e^i
        // c_shares[repetition][party] =
        // eval(P_ei[repetition][party], R_es[repetition]);
        c_shares[repetition][party] = dot_product(
            lagrange_polys_evaluated_at_Re_2m2, P_e_shares[repetition][party]);
      }
    }

    // calculate missing shares
    c[repetition] = proof.P_at_R;
    c_shares[repetition][missing_party] = proof.P_at_R;
    auto a_shares_missing = a_shares.get(repetition, missing_party);
    auto b_shares_missing = b_shares.get(repetition, missing_party);
    for (size_t j = 0; j < instance.m1; j++) {
      a[repetition][j] = proof.S_j_at_R[j];
      a_shares_missing[j] = proof.S_j_at_R[j];
      b[repetition][j] = proof.T_j_at_R[j];
      b_shares_missing[j] = proof.T_j_at_R[j];
    }
    for (size_t party = 0; party < instance.num_MPC_parties; party++) {
      if (party != missing_party) {
        c_shares[repetition][missing_party] -= c_shares[repetition][party];
        auto a_shares_party = a_shares.get(repetition, party);
        auto b_shares_party = b_shares.get(repetition, party);
        for (size_t j = 0; j < instance.m1; j++) {
          a_shares_missing[j] -= a_shares_party[j];
          b_shares_missing[j] -= b_shares_party[j];
        }
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  // recompute h_1 and h_3
  /////////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<uint8_t>> sk_deltas;
  std::vector<std::vector<uint8_t>> t_deltas;
  for (const banquet_repetition_proof_t &proof : signature.proofs) {
    sk_deltas.push_back(proof.sk_delta);
    t_deltas.push_back(proof.t_delta);
  }
  std::vector<uint8_t> h_1 = phase_1_commitment(
      instance, signature.salt, pk, message, message_len,
      party_seed_commitments, sk_deltas, t_deltas, rep_output_broadcasts);

  std::vector<uint8_t> h_3 = phase_3_commitment(
      instance, signature.salt, h_2, c, c_shares, a, a_shares, b, b_shares);
  // do checks
  if (memcmp(h_1.data(), signature.h_1.data(), h_1.size()) != 0) {
    return false;
  }
  if (memcmp(h_2.data(), signature.h_2.data(), h_2.size()) != 0) {
    return false;
  }
  if (memcmp(h_3.data(), signature.h_3.data(), h_3.size()) != 0) {
    return false;
  }

  // check if P_e(R) = Sum_j S_e_j(R) * T_e_j(R) for all e
  for (size_t repetition = 0; repetition < instance.num_rounds; repetition++) {
    field::GF2E accum;
    for (size_t j = 0; j < instance.m1; j++) {
      accum += signature.proofs[repetition].S_j_at_R[j] *
               signature.proofs[repetition].T_j_at_R[j];
    }
    if (accum != signature.proofs[repetition].P_at_R) {
      return false;
    }
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
    for (const std::vector<uint8_t> &seed : proof.reveallist.first) {
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
      std::vector<uint8_t> buffer = proof.P_delta[k].to_bytes();
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
    {
      std::vector<uint8_t> buffer = proof.P_at_R.to_bytes();
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer = proof.S_j_at_R[j].to_bytes();
      serialized.insert(serialized.end(), buffer.begin(), buffer.end());
    }
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer = proof.T_j_at_R[j].to_bytes();
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
  std::vector<uint8_t> h_1(instance.digest_size), h_2(instance.digest_size),
      h_3(instance.digest_size);
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
      std::vector<uint8_t> seed(instance.seed_size);
      memcpy(seed.data(), serialized.data() + current_offset, seed.size());
      current_offset += seed.size();
      reveallist.first.push_back(seed);
    }
    std::vector<uint8_t> C_e(instance.digest_size);
    memcpy(C_e.data(), serialized.data() + current_offset, C_e.size());
    current_offset += C_e.size();

    std::vector<uint8_t> sk_delta(instance.aes_params.key_size);
    memcpy(sk_delta.data(), serialized.data() + current_offset,
           sk_delta.size());
    current_offset += sk_delta.size();

    std::vector<uint8_t> t_delta(instance.aes_params.num_sboxes);
    memcpy(t_delta.data(), serialized.data() + current_offset, t_delta.size());
    current_offset += t_delta.size();

    std::vector<uint8_t> output_broadcast(instance.aes_params.block_size *
                                          instance.aes_params.num_blocks);
    memcpy(output_broadcast.data(), serialized.data() + current_offset,
           output_broadcast.size());
    current_offset += output_broadcast.size();

    field::GF2E tmp;
    std::vector<field::GF2E> P_delta;
    P_delta.reserve(instance.m2 + 1);
    for (size_t k = 0; k < instance.m2 + 1; k++) {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      tmp.from_bytes(buffer.data());
      P_delta.push_back(tmp);
    }
    field::GF2E P_at_R;
    {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      P_at_R.from_bytes(buffer.data());
    }
    std::vector<field::GF2E> S_j_at_R;
    S_j_at_R.reserve(instance.m1);
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      tmp.from_bytes(buffer.data());
      S_j_at_R.push_back(tmp);
    }
    std::vector<field::GF2E> T_j_at_R;
    T_j_at_R.reserve(instance.m1);
    for (size_t j = 0; j < instance.m1; j++) {
      std::vector<uint8_t> buffer(instance.lambda);
      memcpy(buffer.data(), serialized.data() + current_offset, buffer.size());
      current_offset += buffer.size();
      tmp.from_bytes(buffer.data());
      T_j_at_R.push_back(tmp);
    }
    proofs.emplace_back(banquet_repetition_proof_t{
        reveallist, C_e, sk_delta, t_delta, output_broadcast, P_delta, P_at_R,
        S_j_at_R, T_j_at_R});
  }
  assert(current_offset == serialized.size());
  banquet_signature_t signature{salt, h_1, h_2, h_3, proofs};
  return signature;
}
