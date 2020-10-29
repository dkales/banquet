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

#include <cstdint>
#include <cstdlib>

/** Parameter set names */
enum banquet_params_t {
  PARAMETER_SET_INVALID = 0,
  Banquet_L1_Param1 = 1,
  PARAMETER_SET_MAX_INDEX = 2
};

struct banquet_instance_t {

  uint32_t digest_size;     /* bytes */
  uint32_t seed_size;       /* bytes */
  uint32_t num_rounds;      // T
  uint32_t num_MPC_parties; // N

  uint32_t input_size;  /* bytes */
  uint32_t output_size; /* bytes */
  uint32_t view_size;   /* bytes */

  banquet_params_t params;
};

const banquet_instance_t *banquet_instance_get(banquet_params_t param);

#endif
