#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-mno-avx"

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
    export EXTRAFLAGS=
else
    export CC=/usr/bin/clang-5.0
    export CXX=/usr/bin/clang++-5.0
    # Boost 1.55 is too old, override Boost.PP detection of variadic macro support
    export EXTRAFLAGS="-DBOOST_PP_VARIADICS=1"
fi

echo $($CC --version)
echo $($CXX --version)

# list the prebuilt prereqs
ls ${INSTALL_PREFIX}

# where to install MPQC4 (need for testing installed code)
export INSTALL_DIR=${INSTALL_PREFIX}/mpqc4

# make build dir
cd ${BUILD_PREFIX}
mkdir -p mpqc4
cd mpqc4

cmake ${TRAVIS_BUILD_DIR} \
    -DTiledArray_DIR="${INSTALL_PREFIX}/TA/lib/cmake/tiledarray" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_CXX_FLAGS="-ftemplate-depth=1024 -Wno-unused-command-line-argument ${EXTRAFLAGS}" \
    -DLIBINT2_INSTALL_DIR="${INSTALL_PREFIX}/libint" \
    -DMPQC_VALIDATION_TEST_PRINT=true

### build
make -j1 mpqc
### test within build tree
# prepend "setarch `uname -m` -R" if switched back to container builds
make -j1 check
### install and test dev samples
make install
cd ${INSTALL_DIR}/share/doc/mpqc*/examples
cd mp2
  cmake . -DCMAKE_BUILD_TYPE=$BUILD_TYPE
  make mp2
  # prepend "setarch `uname -m` -R" if switched back to container builds
  ./mp2 ./mp2.json
cd ..
