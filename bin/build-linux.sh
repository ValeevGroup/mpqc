#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-mno-avx"

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    export CC=/usr/bin/clang-5.0
    export CXX=/usr/bin/clang++-5.0
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
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/mpqc4" \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_CXX_FLAGS="-ftemplate-depth=1024 -Wno-unused-command-line-argument" \
    -DLIBINT2_INSTALL_DIR="$INSTALL_DIR/libint" \
    -DMPQC_VALIDATION_TEST_PRINT=true

### build
make -j1 mpqc
### test within build tree
setarch `uname -m` -R make -j1 check
### install and test dev samples
make install
cd $INSTALL_DIR/mpqc4/share/doc/mpqc*/examples
cd mp2
  cmake .
  make mp2
  setarch `uname -m` -R ./mp2 ./mp2.json
cd ..
