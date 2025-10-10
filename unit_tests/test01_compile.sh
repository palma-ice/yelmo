!#/bin/bash

# Compile with debugging, for now yes
debug=1

# Clean all compilations
make clean

# Compile static library with debugging on
make yelmo-static debug=$debug

