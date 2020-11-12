/*
 *  This file is part of the optimized implementation of the Picnic signature
 * scheme. See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../banquet.h"
#include "bench_timing.h"
#include "bench_utils.h"

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>

struct timing_and_size_t {
  uint64_t keygen, sign, serialize, deserialize, verify, size;
};

static void print_timings(const std::vector<timing_and_size_t> &timings) {
  printf("keygen,sign,verify,size,serialize,deserialize\n");
  for (const auto &timing : timings) {
    printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64
           ",%" PRIu64 "\n",
           timing.keygen, timing.sign, timing.verify, timing.size,
           timing.serialize, timing.deserialize);
  }
}

static void bench_sign_and_verify_free(const bench_options_free_t *options) {
  static const uint8_t m[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                              12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                              23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

  std::vector<timing_and_size_t> timings(options->iter);

  timing_context_t ctx;
  if (!timing_init(&ctx)) {
    printf("Failed to initialize timing functionality.\n");
    return;
  }
  banquet_instance_t instance;
  instance.digest_size = 2 * options->kappa;
  instance.seed_size = options->kappa;
  switch (options->kappa) {
  case 16:
    instance.aes_params = {16, 16, 1, 200};
    break;
  case 24:
    instance.aes_params = {24, 16, 2, 416};
    break;
  case 32:
    instance.aes_params = {32, 16, 2, 500};
    break;
  default:
    printf("invalid kappa, choose 16,24,32\n");
  }
  instance.m1 = options->m1;
  instance.m2 = options->m2;
  instance.lambda = options->lambda;
  instance.num_MPC_parties = options->N;
  instance.num_rounds = options->tau;
  printf(
      "Instance: N=%d, tau=%d, lambda=%d, m1=%d, m2=%d, AES-Keylen=Seclvl=%d\n",
      instance.num_MPC_parties, instance.num_rounds, instance.lambda,
      instance.m1, instance.m2, instance.aes_params.key_size);
  if (instance.m1 * instance.m2 != instance.aes_params.num_sboxes) {
    printf("invalid m1 and m2 for chosen seclvl\n");
    return;
  }

  for (unsigned int i = 0; i != options->iter; ++i) {
    timing_and_size_t &timing = timings[i];

    uint64_t start_time = timing_read(&ctx);
    banquet_keypair_t keypair = banquet_keygen(instance);

    uint64_t tmp_time = timing_read(&ctx);
    timing.keygen = tmp_time - start_time;
    start_time = timing_read(&ctx);

    banquet_signature_t signature =
        banquet_sign(instance, keypair, m, sizeof(m));

    tmp_time = timing_read(&ctx);
    timing.sign = tmp_time - start_time;
    start_time = timing_read(&ctx);
    std::vector<uint8_t> serialized =
        banquet_serialize_signature(instance, signature);
    tmp_time = timing_read(&ctx);
    timing.serialize = tmp_time - start_time;
    timing.size = serialized.size();

    start_time = timing_read(&ctx);
    banquet_signature_t deserialized =
        banquet_deserialize_signature(instance, serialized);
    tmp_time = timing_read(&ctx);
    timing.deserialize = tmp_time - start_time;
    start_time = timing_read(&ctx);
    bool ok =
        banquet_verify(instance, keypair.second, deserialized, m, sizeof(m));
    tmp_time = timing_read(&ctx);
    timing.verify = tmp_time - start_time;
    if (!ok)
      std::cerr << "failed to verify signature" << std::endl;
  }

  timing_close(&ctx);
  print_timings(timings);
}

int main(int argc, char **argv) {
  bench_options_free_t opts = {0, 0, 0, 0, 0, 0, 0};
  int ret = parse_args_free(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_sign_and_verify_free(&opts);
  }

  return ret;
}
