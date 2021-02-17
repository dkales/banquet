#pragma once

#include "banquet.h"
extern "C" {
#include "kdf_shake.h"
}
#include "gsl-lite.hpp"
#include "tree.h"
#include <cstdlib>

class RandomTape {
private:
  /* data */
  hash_context ctx;

public:
  RandomTape(const gsl::span<uint8_t> &seed, const banquet_salt_t &salt,
             size_t rep_index, size_t party_index);
  ~RandomTape() = default;

  void squeeze_bytes(uint8_t *out, size_t len);
};

class RandomTapes {
private:
  /* data */
  RepByteContainer random_tapes;
  size_t random_tape_size;

public:
  RandomTapes(size_t num_repetitions, size_t num_parties,
              size_t random_tape_size)
      : random_tapes(num_repetitions, num_parties, random_tape_size),
        random_tape_size(random_tape_size){};
  ~RandomTapes() = default;

  void generate_4_tapes(size_t repetition, size_t start_party,
                        const banquet_salt_t &salt,
                        const gsl::span<uint8_t> &seed0,
                        const gsl::span<uint8_t> &seed1,
                        const gsl::span<uint8_t> &seed2,
                        const gsl::span<uint8_t> &seed3);
  void generate_tape(size_t repetition, size_t party,
                     const banquet_salt_t &salt,
                     const gsl::span<uint8_t> &seed);
  gsl::span<uint8_t> get_bytes(size_t repetition, size_t party, size_t start,
                               size_t len);
};