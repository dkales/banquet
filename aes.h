#ifndef _AES_H_
#define _AES_H_

#include "picnic2_types.h"
#include <stdint.h>

int aes_128(const mzd_local_t* key, const mzd_local_t* plaintext, mzd_local_t* ciphertext);
int aes_128_old(const uint8_t* key, const uint8_t* plaintext, uint8_t* ciphertext);
void aux_aes_128(const mzd_local_t* key, randomTape_t* tapes);
int mpc_aes_128(mzd_local_t* maskedKey, shares_aes_t* mask_shares, randomTape_t* tapes, msgs_t* msgs,
                const uint8_t* plaintext, const uint8_t* pubKey, const picnic_instance_t* params);

#endif
