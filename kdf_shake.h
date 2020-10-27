/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef KDF_SHAKE_H
#define KDF_SHAKE_H

#include <stdint.h>

#include "macros.h"
#include "endian_compat.h"

#include "keccak/KeccakHash.h"

typedef Keccak_HashInstance hash_context ATTR_ALIGNED(32);

/**
 * Initialize hash context based on the digest size used by Picnic. If the size is 32 bytes,
 * SHAKE128 is used, otherwise SHAKE256 is used.
 */
static inline void hash_init(hash_context *ctx, size_t digest_size)
{
  if (digest_size == 32)
  {
    Keccak_HashInitialize_SHAKE128(ctx);
  }
  else
  {
    Keccak_HashInitialize_SHAKE256(ctx);
  }
}

static inline void hash_update(hash_context *ctx, const uint8_t *data, size_t size)
{
  Keccak_HashUpdate(ctx, data, size << 3);
}

static inline void hash_final(hash_context *ctx)
{
  Keccak_HashFinal(ctx, NULL);
}

static inline void hash_squeeze(hash_context *ctx, uint8_t *buffer, size_t buflen)
{
  Keccak_HashSqueeze(ctx, buffer, buflen << 3);
}

#define hash_clear(ctx)

static inline void hash_update_uint16_le(hash_context *ctx, uint16_t data)
{
  const uint16_t data_le = htole16(data);
  hash_update(ctx, (const uint8_t *)&data_le, sizeof(data_le));
}

static inline void hash_init_prefix(hash_context *ctx, size_t digest_size,
                                    const uint8_t prefix)
{
  hash_init(ctx, digest_size);
  hash_update(ctx, &prefix, sizeof(prefix));
}

typedef hash_context kdf_shake_t;

#define kdf_shake_init(ctx, digest_size) hash_init((ctx), (digest_size))
#define kdf_shake_init_prefix(ctx, digest_size, prefix) hash_init_prefix((ctx), (digest_size), (prefix))
#define kdf_shake_update_key(ctx, key, keylen) hash_update((ctx), (key), (keylen))
#define kdf_shake_update_key_uint16_le(ctx, key) hash_update_uint16_le((ctx), (key))
#define kdf_shake_finalize_key(ctx) hash_final((ctx))
#define kdf_shake_get_randomness(ctx, dst, count) hash_squeeze((ctx), (dst), (count))
#define kdf_shake_clear(ctx) hash_clear((ctx))

#endif