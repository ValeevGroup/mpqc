#! /bin/sh

# set to the release id of the required library
export RELID=2.3.0-beta.2

# Exit on error
set -ev
case "$CXX" in
    g++)
        export CC=/usr/bin/gcc-$GCC_VERSION
        export CXX=/usr/bin/g++-$GCC_VERSION
        ;;
    clang++)
        export CC=/usr/bin/clang-3.8
        export CXX=/usr/bin/clang++-3.8
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

ls 
pwd

# Make install directory for MPQC dependencies 
mkdir -p /home/travis/build/ValeevGroup/_build
mkdir -p /home/travis/build/ValeevGroup/_install

cd /home/travis/build/ValeevGroup/_build

# download and unpack libint tarball
wget --no-check-certificate -q https://github.com/evaleev/libint/releases/download/v$RELID/libint-$RELID-test-mpqc4.tgz
tar -xvzf libint-$RELID-test-mpqc4.tgz
cd libint-$RELID/

./configure --prefix="/home/travis/build/ValeevGroup/_install/libint" \
 --with-incdirs="-I/usr/include/eigen3"


make -j2
make install
