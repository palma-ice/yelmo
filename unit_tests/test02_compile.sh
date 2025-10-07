!#/bin/bash

# This test should be run after test01_compile, or run `make clean` beforehand

# Run with debugging, for now yes
debug=1

# Make available programs
make benchmarks debug=$debug

make initmip debug=$debug

make trough debug=$debug

make ismiphom debug=$debug

make slab debug=$debug

make mismip debug=$debug
