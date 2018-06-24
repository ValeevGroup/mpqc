#! /bin/sh

# set to the release id of the required library
export RELID=2.4.1

# Exit on error
set -ev
case "$CXX" in
    g++)
        export CC=/usr/bin/gcc-$GCC_VERSION
        export CXX=/usr/bin/g++-$GCC_VERSION
        ;;
    clang++)
        export CC=/usr/bin/clang-$CLANG_VERSION
        export CXX=/usr/bin/clang++-$CLANG_VERSION
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

# re-build if not in cache already ... must manualy wipe cache when need to bump version or update toolchain
export INSTALL_DIR=${INSTALL_PREFIX}/libint
if [ ! -d "${INSTALL_DIR}" ]; then
  cd ${BUILD_PREFIX}

  # download and unpack libint tarball
  wget --no-check-certificate -q https://github.com/evaleev/libint/releases/download/v$RELID/libint-$RELID-test-mpqc4.tgz
  tar -xvzf libint-$RELID-test-mpqc4.tgz
  cd libint-$RELID/

  ./configure --prefix="${INSTALL_DIR}" \
    --with-incdirs="-I/usr/include/eigen3" --disable-shared --enable-static

  make -j2
  make install
fi
