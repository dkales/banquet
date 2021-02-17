#include "tape.h"

RandomTape::RandomTape(const gsl::span<uint8_t> &seed,
                       const banquet_salt_t &salt, size_t rep_index,
                       size_t party_index) {
  hash_init(&ctx, seed.size() * 2);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)rep_index);
  hash_update_uint16_le(&ctx, (uint16_t)party_index);
  hash_final(&ctx);
}

void RandomTape::squeeze_bytes(uint8_t *out, size_t len) {
  hash_squeeze(&ctx, out, len);
}

void RandomTapes::generate_4_tapes(size_t repetition, size_t start_party,
                                   const banquet_salt_t &salt,
                                   const gsl::span<uint8_t> &seed0,
                                   const gsl::span<uint8_t> &seed1,
                                   const gsl::span<uint8_t> &seed2,
                                   const gsl::span<uint8_t> &seed3) {

  hash_context_x4 ctx;
  hash_init_x4(&ctx, seed0.size() * 2);
  hash_update_x4_4(&ctx, seed0.data(), seed1.data(), seed2.data(), seed3.data(),
                   seed0.size());
  hash_update_x4_1(&ctx, salt.data(), salt.size());
  hash_update_x4_uint16_le(&ctx, (uint16_t)repetition);
  const uint16_t parties[4] = {
      (uint16_t)(start_party), (uint16_t)(start_party + 1),
      (uint16_t)(start_party + 2), (uint16_t)(start_party + 3)};
  hash_update_x4_uint16s_le(&ctx, parties);
  hash_final_x4(&ctx);
  hash_squeeze_x4_4(&ctx, random_tapes.get(repetition, parties[0]).data(),
                    random_tapes.get(repetition, parties[1]).data(),
                    random_tapes.get(repetition, parties[2]).data(),
                    random_tapes.get(repetition, parties[3]).data(),
                    random_tape_size);
}
void RandomTapes::generate_tape(size_t repetition, size_t party,
                                const banquet_salt_t &salt,
                                const gsl::span<uint8_t> &seed) {
  hash_context ctx;
  hash_init(&ctx, seed.size() * 2);
  hash_update(&ctx, seed.data(), seed.size());
  hash_update(&ctx, salt.data(), salt.size());
  hash_update_uint16_le(&ctx, (uint16_t)repetition);
  hash_update_uint16_le(&ctx, (uint16_t)party);
  hash_final(&ctx);
  hash_squeeze(&ctx, random_tapes.get(repetition, party).data(),
               random_tape_size);
}

gsl::span<uint8_t> RandomTapes::get_bytes(size_t repetition, size_t party,
                                          size_t start, size_t len) {
  auto tape = random_tapes.get(repetition, party);
  return tape.subspan(start, len);
}