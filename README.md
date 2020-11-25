# banquet
Banquet implementation



## Requirements

* [GMP](https://gmplib.org/)
* [NTL](https://shoup.net/ntl)

It is recommended to build NTL yourself since the shared library version shipped by some distros have a 2x performance difference.

## Setup

```bash
mkdir build
cd build
cmake ..
make 
# tests (if you built them)
make test
# benchmarks
./bench -i <iterations> <instance> #instance from banquet_instances.h
./bench_free -i <iterations> <kappa> <N> <tau> <m1> <m2> <lambda> #benchmark parameters freely
```
