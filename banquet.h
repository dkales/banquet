#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "banquet_instances.h"
#include "types.h"

// crypto api
banquet_keypair_t banquet_keygen(const banquet_instance_t &instance);

banquet_signature_t banquet_sign(const banquet_instance_t &instance,
                                 const banquet_keypair_t &keypair,
                                 const uint8_t *message, size_t message_len);