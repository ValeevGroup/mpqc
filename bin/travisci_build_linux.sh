#!/bin/sh

set -e

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
fi
export FC=/usr/bin/gfortran-$GCC_VERSION
export CXXFLAGS="-std=c++11"
export LDFLAGS=$OPENMPFLAGS

./configure
make -j2
make check0
