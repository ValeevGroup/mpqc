#! /bin/sh

${TRAVIS_BUILD_DIR}/bin/build-mpich-$TRAVIS_OS_NAME.sh
${TRAVIS_BUILD_DIR}/bin/build-libint-$TRAVIS_OS_NAME.sh
${TRAVIS_BUILD_DIR}/bin/build-TA-$TRAVIS_OS_NAME.sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-mno-avx"
echo $GCC_VERSION
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
    export EXTRAFLAGS="-DHAVE_LAPACK_CONFIG_H=1 -DLAPACK_COMPLEX_CPP=1 -DBTAS_HAS_CBLAS=1"
else
    export CC=/usr/bin/clang-6.0
    export CXX=/usr/bin/clang++-6.0
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


if [ "$BUILD_TYPE" = "Debug" ] && [ "$GCC_VERSION" = 5 ]; then
    export CODECOVCXXFLAGS="--coverage -O0"
fi

cmake ${TRAVIS_BUILD_DIR} \
    -DTiledArray_DIR="${INSTALL_PREFIX}/TA/lib/cmake/tiledarray" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_CXX_FLAGS="-ftemplate-depth=1024 -Wno-unused-command-line-argument ${EXTRAFLAGS} ${CODECOVCXXFLAGS}" \
    -DLIBINT2_INSTALL_DIR="${INSTALL_PREFIX}/libint" \
    -DMPQC_VALIDATION_TEST_PRINT=true \
    -DCMAKE_TOOLCHAIN_FILE="/home/travis/build/ValeevGroup/mpqc4/bin/travis-lapacke.cmake"

### build
make -j1 mpqc VERBOSE=1
### test within build tree
setarch `uname -m` -R make -j1 check
### install and test dev samples
make install
cd ${INSTALL_DIR}/share/doc/mpqc*/examples
cd mp2
  cmake . \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DCMAKE_CXX_FLAGS=$CODECOVCXXFLAGS
  make mp2
  setarch `uname -m` -R ./mp2 ./mp2.json
cd ..

#    -Dlapack_complex_float="std::complex<float>" \
#    -Dlapack_complex_double="std::complex<double>" \
