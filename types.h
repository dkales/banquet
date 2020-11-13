#pragma once

#include "gsl-lite.hpp"
#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "field.h"

/* Prefix values for domain separation */
constexpr uint8_t HASH_PREFIX_0 = 0;
constexpr uint8_t HASH_PREFIX_1 = 1;
constexpr uint8_t HASH_PREFIX_2 = 2;
constexpr uint8_t HASH_PREFIX_3 = 3;
constexpr uint8_t HASH_PREFIX_4 = 4;
constexpr uint8_t HASH_PREFIX_5 = 5;

constexpr size_t SALT_SIZE = 32;
typedef std::array<uint8_t, SALT_SIZE> banquet_salt_t;

typedef std::pair<std::vector<std::vector<uint8_t>>, size_t> reveal_list_t;

typedef std::array<uint8_t, 16> aes_block_t;

typedef std::pair<std::vector<uint8_t>, std::vector<uint8_t>> banquet_keypair_t;

struct banquet_repetition_proof_t {
  reveal_list_t reveallist;
  std::vector<uint8_t> C_e;
  std::vector<uint8_t> sk_delta;
  std::vector<uint8_t> t_delta;
  std::vector<field::GF2E> P_delta;
  field::GF2E P_at_R;
  std::vector<field::GF2E> S_j_at_R;
  std::vector<field::GF2E> T_j_at_R;
};

struct banquet_signature_t {
  banquet_salt_t salt;
  std::vector<uint8_t> h_1;
  std::vector<uint8_t> h_3;
  std::vector<banquet_repetition_proof_t> proofs;
};

template <typename T> class RepContainer {
  std::vector<T> _data;
  size_t _num_repetitions;
  size_t _num_parties;
  size_t _object_size;

public:
  RepContainer(size_t num_repetitions, size_t num_parties, size_t object_size)
      : _data(num_repetitions * num_parties * object_size),
        _num_repetitions(num_repetitions), _num_parties(num_parties),
        _object_size(object_size) {}

  inline gsl::span<T> get(size_t repetition, size_t party) {
    size_t offset =
        (repetition * _num_parties * _object_size) + (party * _object_size);
    return gsl::span<T>(_data.data() + offset, _object_size);
  }
  inline gsl::span<const T> get(size_t repetition, size_t party) const {
    size_t offset =
        (repetition * _num_parties * _object_size) + (party * _object_size);
    return gsl::span<const T>(_data.data() + offset, _object_size);
  }

  std::vector<gsl::span<T>> get_repetition(size_t repetition) {
    std::vector<gsl::span<T>> ret;
    ret.reserve(_num_parties);
    size_t offset = (repetition * _num_parties * _object_size);
    for (size_t i = 0; i < _num_parties; i++)
      ret.emplace_back(_data.data() + offset + i * _object_size, _object_size);
    return ret;
  }
};

typedef RepContainer<uint8_t> RepByteContainer;