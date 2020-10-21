/*
 *  This file is part of the optimized implementation of the BANQUET signature
 * scheme. See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef BANQUET_INSTANCES_H
#define BANQUET_INSTANCES_H

#include "banquet.h"
#include <stdint.h>

#define SALT_SIZE 32
#define MAX_DIGEST_SIZE 64

typedef struct banquet_instance_t {

  uint32_t digest_size;     /* bytes */
  uint32_t seed_size;       /* bytes */
  uint32_t num_rounds;      // T
  uint32_t num_MPC_parties; // N

  uint32_t input_size;  /* bytes */
  uint32_t output_size; /* bytes */
  uint32_t view_size;   /* bytes */

  banquet_params_t params;

} banquet_instance_t;

const banquet_instance_t *banquet_instance_get(banquet_params_t param);

/* Prefix values for domain separation */
static const uint8_t HASH_PREFIX_0 = 0;
static const uint8_t HASH_PREFIX_1 = 1;
static const uint8_t HASH_PREFIX_2 = 2;
static const uint8_t HASH_PREFIX_3 = 3;
static const uint8_t HASH_PREFIX_4 = 4;
static const uint8_t HASH_PREFIX_5 = 5;

#endif
