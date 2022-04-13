#ifndef BANQUET_INSTANCES_H
#define BANQUET_INSTANCES_H

#include <cstdint>
#include <cstdlib>

/** Parameter set names */
enum banquet_params_t {
  PARAMETER_SET_INVALID = 0,
  Banquet_L0_Param1 = 1,
  Banquet_L1_Param1 = 2,
  Banquet_L1_Param2 = 3,
  Banquet_L1_Param3 = 4,
  Banquet_L1_Param4 = 5,
  Banquet_L1_Param5 = 6,
  Banquet_L1_Param6 = 7,
  Banquet_L1_Param7 = 8,
  Banquet_L1_Param8 = 9,
  Banquet_L1_Param9 = 10,
  Banquet_L1_Param10 = 11,
  Banquet_L3_Param1 = 12,
  Banquet_L3_Param2 = 13,
  Banquet_L3_Param3 = 14,
  Banquet_L3_Param4 = 15,
  Banquet_L3_Param5 = 16,
  Banquet_L3_Param6 = 17,
  Banquet_L3_Param7 = 18,
  Banquet_L5_Param1 = 19,
  Banquet_L5_Param2 = 10,
  Banquet_L5_Param3 = 21,
  Banquet_L5_Param4 = 22,
  Banquet_L5_Param5 = 23,
  Banquet_L5_Param6 = 24,
  PARAMETER_SET_MAX_INDEX = 25
};

struct banquet_aes_t {
  uint32_t key_size;
  uint32_t block_size;
  uint32_t num_blocks;
  uint32_t num_sboxes;
};

struct banquet_instance_t {

  banquet_aes_t aes_params;

  uint32_t digest_size;     /* bytes */
  uint32_t seed_size;       /* bytes */
  uint32_t num_rounds;      // T
  uint32_t num_MPC_parties; // N

  // m1 * m2 = m (NUM_AES_SBOXES)
  uint32_t m1;     // m1: dimension 1 of sqrt check
  uint32_t m2;     // m2: dimension 2 of sqrt check
  uint32_t lambda; // field expansion factor

  banquet_params_t params;
};

const banquet_instance_t &banquet_instance_get(banquet_params_t param);

#endif
