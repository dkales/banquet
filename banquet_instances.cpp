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

static banquet_instance_t instances[PARAMETER_SET_MAX_INDEX] = {
    {
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        PARAMETER_SET_INVALID,
    },
    /* digest size, seed size, T, N, input size, output size, view size, m1, m2, lambda */
    {32, 16, 31, 64, 16, 16, 0, 10, 20, 4, Banquet_L1_Param1},
    {32, 16, 31, 64, 16, 16, 0, 20, 10, 4, Banquet_L1_Param2},
    {32, 16, 29, 64, 16, 16, 0, 10, 20, 5, Banquet_L1_Param3},
    {32, 16, 27, 64, 16, 16, 0, 10, 20, 6, Banquet_L1_Param4},
    {32, 16, 28, 128, 16, 16, 0, 10, 20, 4, Banquet_L1_Param5},
    {32, 16, 26, 128, 16, 16, 0, 10, 20, 5, Banquet_L1_Param6},
    {32, 16, 24, 128, 16, 16, 0, 10, 20, 6, Banquet_L1_Param7},
    {32, 16, 25, 256, 16, 16, 0, 10, 20, 4, Banquet_L1_Param8},
    {32, 16, 23, 256, 16, 16, 0, 10, 20, 5, Banquet_L1_Param9},
    {32, 16, 21, 256, 16, 16, 0, 10, 20, 6, Banquet_L1_Param10},
};

const banquet_instance_t &banquet_instance_get(banquet_params_t param) {
  if (param <= PARAMETER_SET_INVALID || param >= PARAMETER_SET_MAX_INDEX) {
    throw std::runtime_error("invalid parameter set");
  }

  return instances[param];
}
