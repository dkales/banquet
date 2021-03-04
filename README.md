# Banquet: Short and Fast Signatures from AES 
This an implementation of the Banquet signature scheme presented in the following work ([eprint](https://eprint.iacr.org/2021/068)):

**Banquet: Short and Fast Signatures from AES**
*Carsten Baum, Cyprien Delpech de Saint Guilhem, Daniel Kales, Emmanuela Orsini, Peter Scholl and Greg Zaverucha*,
Public Key Cryptography 2021


Banquet is a signature scheme build on the MPC-in-the-Head paradigm and is intended to provide security against both classical and quantum attackers. The security against quantum attackers stems from the use of symmetric-key cryptography such as block ciphers and hash functions, without relying on certain number-theoretic hardness assumptions that are broken in the presence of quantum attackers. In contrast to similar designs such as [Picnic](https://microsoft.github.io/Picnic/), Banquet relies only on standardized symmetric-key primitives such as AES and SHAKE. While these choices traditionally make for much larger and slower signatures, Banquet signatures can nearly match the signature sizes of Picnic (albeit with slower, but still practical run times) or have speed within a factor of two of Picnic (at the cost of larger signatures).


## Requirements

Dependencies for building
* A C++17 compatible C++ tool chain
* [CMake](https://cmake.org/) 3.12+

Additional requirements for unit tests
* [GMP](https://gmplib.org/)
* [NTL](https://shoup.net/ntl)


## Setup

```bash
mkdir build
cd build
cmake ..
make 
# tests (if you built them by passing -DBUILD_TESTS=On to CMake)
make test
# benchmarks
./bench -i <iterations> <instance> #instance from banquet_instances.h
./bench_free -i <iterations> <kappa> <N> <tau> <m1> <m2> <lambda> #benchmark parameters freely, see paper for secure parameters
```

## Acknowledgements

Some parts of the code were based on the [optimized Picnic implementation](https://github.com/IAIK/Picnic).