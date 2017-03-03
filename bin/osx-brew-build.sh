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
#    - libint
wget https://github.com/evaleev/libint/releases/download/v2.3.0-beta.3/libint-2.3.0-beta.3.tgz && tar -xvzf libint-2.3.0-beta.3.tgz && cd libint-2.3.0-beta.3 && ./configure --prefix=$PREFIX --with-incdirs="-I/usr/local/include/eigen3" --enable-shared --disable-static && make -j2 && make install && make clean && cd .. && rm -rf libint-2.3.0-beta.3
#    - tiledarray
mkdir tiledarray && cd tiledarray && git clone --depth=1 https://github.com/ValeevGroup/tiledarray.git tiledarray_src && cmake tiledarray_src -DCMAKE_INSTALL_PREFIX=$PREFIX -DMPI_CXX_COMPILER=mpicxx -DMPI_C_COMPILER=mpicc -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=tiledarray_src/cmake/toolchains/osx-clang-mpi-accelerate.cmake && make && make install && make clean && cd .. && rm -rf tiledarray

# 3. build and install MPQC4
git clone https://evaleev:f0aee3276b87c17a47d8b18e7c82af7a1cad8842@github.com/ValeevGroup/mpqc4.git && mkdir mpqc4-build && cd mpqc4-build && cmake ../mpqc4 -DCMAKE_INSTALL_PREFIX=$PREFIX -DTiledArray_DIR="$PREFIX/lib/cmake/tiledarray" -DLIBINT2_INSTALL_DIR=$PREFIX -DCMAKE_BUILD_TYPE=Release && make mpqc && make install

##############################################################

##############################################################
# done
clean_up
