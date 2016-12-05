#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-mno-avx -std=c++11"

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    export CC=/usr/bin/clang-3.8
    export CXX=/usr/bin/clang++-3.8
fi

echo $($CC --version)
echo $($CXX --version)

mkdir -p mpqc4_build
cd mpqc4_build

INSTALL_DIR=/home/travis/build/ValeevGroup/_install
mkdir -p /home/travis/build/ValeevGroup/_install
ls $INSTALL_DIR

cmake .. \
    -DTiledArray_DIR="$INSTALL_DIR/TA/lib/cmake/tiledarray" \
    -DCMAKE_PREFIX_PATH="$INSTALL_DIR/mpqc4" \
    -DCMAKE_BUILD_TYPE=DEBUG \
    -DCMAKE_CXX_FLAGS="-ftemplate-depth=1024 -Wno-unused-command-line-argument" \
    -DLIBINT2_INSTALL_DIR="$INSTALL_DIR/libint"

make -j1 mpqc
make -j1 check
