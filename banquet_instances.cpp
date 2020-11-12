/*
 *  This file is part of the optimized implementation of the Picnic signature
 * scheme. See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include "banquet_instances.h"

#include <stdexcept>

/* key_size, block_size, num_blocks, num_sboxes */
constexpr banquet_aes_t AES128_PARAMS = {16, 16, 1, 200 /* 160 + 40 */};
constexpr banquet_aes_t AES192_PARAMS = {24, 16, 2, 416 /* 2*192 + 32 */};
constexpr banquet_aes_t AES256_PARAMS = {32, 16, 2, 500 /* 2*224 + 52 */};

static const banquet_instance_t instances[PARAMETER_SET_MAX_INDEX] = {
    {
        {0, 0, 0, 0},
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        PARAMETER_SET_INVALID,
    },
    /* AES_params, digest size, seed size, T, N, m1, m2, lambda */
    {AES128_PARAMS, 32, 16, 31, 64, 10, 20, 4, Banquet_L1_Param1},
    {AES128_PARAMS, 32, 16, 31, 64, 20, 10, 4, Banquet_L1_Param2},
    {AES128_PARAMS, 32, 16, 29, 64, 10, 20, 5, Banquet_L1_Param3},
    {AES128_PARAMS, 32, 16, 27, 64, 10, 20, 6, Banquet_L1_Param4},
    {AES128_PARAMS, 32, 16, 41, 16, 10, 20, 4, Banquet_L1_Param5},
    {AES128_PARAMS, 32, 16, 35, 32, 10, 20, 4, Banquet_L1_Param6},
    {AES128_PARAMS, 32, 16, 28, 128, 10, 20, 4, Banquet_L1_Param8},
    {AES128_PARAMS, 32, 16, 24, 128, 10, 20, 6, Banquet_L1_Param7},
    {AES128_PARAMS, 32, 16, 23, 256, 10, 20, 5, Banquet_L1_Param9},
    {AES128_PARAMS, 32, 16, 21, 256, 10, 20, 6, Banquet_L1_Param10},
    {AES192_PARAMS, 48, 24, 46, 64, 16, 26, 4, Banquet_L3_Param1},
    {AES192_PARAMS, 48, 24, 46, 64, 26, 16, 4, Banquet_L3_Param2},
    {AES192_PARAMS, 48, 24, 62, 16, 26, 16, 4, Banquet_L3_Param3},
    {AES192_PARAMS, 48, 24, 53, 32, 26, 16, 4, Banquet_L3_Param4},
    {AES192_PARAMS, 48, 24, 40, 64, 26, 16, 6, Banquet_L3_Param5},
    {AES192_PARAMS, 48, 24, 36, 128, 26, 16, 6, Banquet_L3_Param6},
    {AES192_PARAMS, 48, 24, 32, 256, 26, 16, 6, Banquet_L3_Param7},
    {AES256_PARAMS, 64, 32, 63, 64, 20, 25, 4, Banquet_L5_Param1},
    {AES256_PARAMS, 64, 32, 84, 16, 25, 20, 4, Banquet_L5_Param2},
    {AES256_PARAMS, 64, 32, 63, 32, 25, 20, 6, Banquet_L5_Param3},
    {AES256_PARAMS, 64, 32, 54, 64, 25, 20, 6, Banquet_L5_Param4},
    {AES256_PARAMS, 64, 32, 48, 128, 25, 20, 6, Banquet_L5_Param5},
    {AES256_PARAMS, 64, 32, 43, 256, 25, 20, 6, Banquet_L5_Param6},
};

const banquet_instance_t &banquet_instance_get(banquet_params_t param) {
  if (param <= PARAMETER_SET_INVALID || param >= PARAMETER_SET_MAX_INDEX) {
    throw std::runtime_error("invalid parameter set");
  }

  return instances[param];
}
