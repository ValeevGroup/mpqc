#! /bin/sh

# Exit on error
set -ev
case "$CXX" in
    g++)
        export CC=/usr/bin/gcc-$GCC_VERSION
        export CXX=/usr/bin/g++-$GCC_VERSION
        ;;
    clang++)
        export CC=/usr/bin/clang-3.7
        export CXX=/usr/bin/clang++-3.7
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

export CXXFLAGS="-O0 -std=c++11"

# Make install directory for MPQC dependencies 
mkdir -p _build
mkdir -p _install

cd _build

# Unpack libint tarball
tar -xvzf ./external/libint-2.2.0-beta1.tgz

cd libint-2.2.0-beta1
./autogen.sh

cd ..
mkdir libint

cd libint

../libint-2.2.0-beta1/configure \
    --prefix="$HOME/_install/libint" \
    CXX=$CXX \
    CXXFLAGS=$CXXFLAGS \
    --with-incdirs=/usr/include/eigen3 

make -j2
make install

# Install MPICH
# if [ ! -d "${HOME}/mpich" ]; then
#     wget --no-check-certificate -q http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
#     tar -xzf mpich-3.2.tar.gz
#     cd mpich-3.2
#     ./configure CC=$CC CXX=$CXX --disable-fortran --disable-romio --prefix=${HOME}/mpich
#     make -j2
#     make install
#     ${HOME}/mpich/bin/mpichversion
#     ${HOME}/mpich/bin/mpicc -show
#     ${HOME}/mpich/bin/mpicxx -show
# else
#     echo "MPICH installed..."
#     find ${HOME}/mpich -name mpiexec
#     find ${HOME}/mpich -name mpicc
#     find ${HOME}/mpich -name mpicxx
# fi
