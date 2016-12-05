#! /bin/sh

# Exit on error
set -ev

# Install packages

case "$CXX" in
    g++)
        export CC=/usr/bin/gcc-$GCC_VERSION
        export CXX=/usr/bin/g++-$GCC_VERSION
        ;;
    clang++)
        export CC=/usr/bin/clang-3.8
        export CXX=/usr/bin/clang++-3.8
        export CXXFLAGS="-std=c++11"
        ;;
    *)
        echo "Unknown C++ compiler:"
        echo "$CXX"
        exit 1
        ;;
esac

# Print compiler information
$CC --version
$CXX --version

# log the CMake version (need 3+)
cmake --version

# Install MPICH
export PREFIX=${HOME}/ValeevGroup/_install/mpich
if [ ! -d "${PREFIX}" ]; then
    wget --no-check-certificate -q http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar -xzf mpich-3.2.tar.gz
    cd mpich-3.2
    ./configure CC=$CC CXX=$CXX --disable-fortran --disable-romio --prefix=${PREFIX}
    make -j2
    make install
    ${PREFIX}/bin/mpichversion
    ${PREFIX}/bin/mpicc -show
    ${PREFIX}/bin/mpicxx -show
else
    echo "MPICH installed..."
    find ${PREFIX} -name mpiexec
    find ${PREFIX} -name mpicc
    find ${PREFIX} -name mpicxx
fi
