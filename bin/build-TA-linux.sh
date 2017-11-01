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

export MPI_HOME=${INSTALL_PREFIX}/mpich
export MPICC=$MPI_HOME/bin/mpicc
export MPICXX=$MPI_HOME/bin/mpicxx
export LD_LIBRARY_PATH=/usr/lib/lapack:/usr/lib/libblas:$LD_LIBRARY_PATH

# Install TA unless previous install is cached ... must manually wipe cache on version bump or toolchain update
export INSTALL_DIR=${INSTALL_PREFIX}/TA
if [ ! -d "${INSTALL_DIR}" ]; then

  # Configure TiledArray
  cd ${BUILD_PREFIX}
  mkdir -p TA
  cd TA

  git clone https://github.com/ValeevGroup/tiledarray.git ta_src

  cmake ta_src \
      -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_C_COMPILER=$CC \
      -DMPI_CXX_COMPILER=$MPICXX \
      -DMPI_C_COMPILER=$MPICC \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE

  # Build all libraries, examples, and applications
  make -j2 VERBOSE=1
  make install
fi
