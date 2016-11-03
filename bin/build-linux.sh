#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-mno-avx -std=c++11"

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    export CC=/usr/bin/clang-3.7
    export CXX=/usr/bin/clang++-3.7
fi

echo $($CC --version)
echo $($CXX --version)

cd _build

mkdir -p mpqc4
cd mpqc4

INSTALL_DIR=/home/travis/build/ValeevGroup/mpqc4/_install
ls $INSTALL_DIR

cmake ../.. \
    -DTiledArray_DIR="$INSTALL_DIR/TA/lib/cmake/tiledarray" \
    -DCMAKE_PREFIX_PATH="$INSTALL_DIR/TA" \
    -DCMAKE_BUILD_TYPE=DEBUG \
    -DLIBINT2_INSTALL_DIR="$INSTALL_DIR/libint"

make -j2 mpqc_info

