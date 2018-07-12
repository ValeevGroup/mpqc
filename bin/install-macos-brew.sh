#!/bin/sh

# builds MPQC and all prerequisites on OS X
# uses HomeBrew to pull the basic toolchain

set -e

function clean_up {
  echo "cleaning up"
  exit
}

function check_program {
  command -v $1 >/dev/null 2>&1 || brew install $1
  echo "dependency \"$1\" satisfied"
}

function check_package {
  brew info $1 | grep -q "Not installed" && brew install $1 $2
  echo "dependency \"$1\" satisfied"
}

trap clean_up SIGHUP SIGINT SIGTERM

##############################################################
########## Prerequisites #####################################
# 0. must have homebrew
command -v brew >/dev/null 2>&1 || { echo >&2 "install HomeBrew first (see http://brew.sh) and put it in your PATH"; exit 1; }

# 1. install pre-prerequisites
check_program git
check_program wget
check_program cmake
check_package openmpi "--c++11 --with-mpi-thread-multiple"
check_package boost
check_package eigen
check_package tbb

# 2. build and install prerequisites

# by default: build on 4 procs
if test -z $NPROC; then
  export NPROC=4
fi

# by default: install under ./install
if test -z $PREFIX; then
  export PREFIX=`pwd`/install
fi

# by default: build TA and MPQC4 in debug mode (Libint always in Release)
if test -z $BUILD_TYPE; then
  export BUILD_TYPE=Debug
fi

#    - libint
export LIBINT_RELID=2.5.0-beta.1
if test ! -d build/libint-${LIBINT_RELID}; then
  mkdir -p build
  cd build
  wget https://github.com/evaleev/libint/releases/download/v${LIBINT_RELID}/libint-${LIBINT_RELID}.tgz && tar -xvzf libint-${LIBINT_RELID}.tgz && cd libint-${LIBINT_RELID}
  ./configure --prefix=$PREFIX --with-incdirs="-I/usr/local/include/eigen3" --enable-shared --disable-static
  make -j${NPROC} && make install && make distclean
  cd ../..
fi

#    - tiledarray
if test ! -d tiledarray; then
  git clone --depth=1 https://github.com/ValeevGroup/tiledarray.git tiledarray
fi
if test ! -d build/tiledarray-clang; then
  mkdir -p build/tiledarray-clang && cd build/tiledarray-clang
  cmake ../../tiledarray -DCMAKE_INSTALL_PREFIX=$PREFIX -DMPI_CXX_COMPILER=mpicxx -DMPI_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_TOOLCHAIN_FILE=../../tiledarray/cmake/toolchains/osx-clang-mpi-accelerate.cmake
  make -j${NPROC} && make install
  cd ../..
fi

# 3. build and install MPQC4
if test ! -d mpqc4; then
  git clone https://evaleev:f0aee3276b87c17a47d8b18e7c82af7a1cad8842@github.com/ValeevGroup/mpqc4.git
fi
if test ! -d build/mpqc4-clang; then
  mkdir -p build/mpqc4-clang && cd build/mpqc4-clang
  cmake ../../mpqc4 -DCMAKE_INSTALL_PREFIX=$PREFIX -DTiledArray_DIR="$PREFIX/lib/cmake/tiledarray" -DLIBINT2_INSTALL_DIR=$PREFIX -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DCMAKE_TOOLCHAIN_FILE=../../tiledarray/cmake/toolchains/osx-clang-mpi-accelerate.cmake
  make -j${NPROC} mpqc && make install
  cd ../..
fi

##############################################################

##############################################################
# done
clean_up
