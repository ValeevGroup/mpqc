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

ls 

# Make install directory for MPQC dependencies 
mkdir -p _build
mkdir -p _install

cd _build

# Unpack libint tarball
tar -xvzf ../external/libint-2.2.0-beta1.tgz
cd libint-2.2.0-beta1/

ls /usr/include
ls /usr/local/include

 ./configure \
     --prefix="$HOME/_install/libint" \
     --with-incdirs=/usr/include/eigen3 \
     --with-cxx-optflags="-O0" 


make -j2
make install

cd ../..
ls 
ls _install
