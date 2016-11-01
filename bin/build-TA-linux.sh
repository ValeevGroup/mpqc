#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-mno-avx"

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    export CC=/usr/bin/clang-3.7
    export CXX=/usr/bin/clang++-3.7
fi

echo $($CC --version)
echo $($CXX --version)

export MPICC=$HOME/mpich/bin/mpicc
export MPICXX=$HOME/mpich/bin/mpicxx
export LD_LIBRARY_PATH=/usr/lib/lapack:/usr/lib/libblas:$LD_LIBRARY_PATH

# Configure TiledArray

mkdir -p _build
cd _build

mkdir -p TA
cd TA

git clone https://github.com/ValeevGroup/tiledarray.git ta_src

cmake ta_src \
      -DCMAKE_INSTALL_PREFIX=../../_install/TA \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_C_COMPILER=$CC \
      -DMPI_CXX_COMPILER=$MPICXX \
      -DMPI_C_COMPILER=$MPICC \
      -DCMAKE_BUILD_TYPE=Debug 

# Build all libraries, examples, and applications
make -j2 VERBOSE=1
make install
