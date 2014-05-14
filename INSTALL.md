To compile MPQC using cmake:

CMAKE 2.8.8 or above is required.
If no cmake is found in path, one will be downloaded and built in ./cmake
when first running ./configure.

Running
$ ./configure cmake 
will FORCE to build and use CMAKE in ./cmake

MPQC relies on several external packages. Some of them we build in ./external directory:
* Eigen is bundled in ./external/eigen.tar.gz; it will be unpacked if needed.
* Boost is bundled in ./external/boost*.tar.bz2; it will be built regardless of whether
  there is one is installed, since Boost depends heavily on compiler.
* LAPACK and HDF5 will be downloaded and built only if none are found on the system.
* Psi3 is bundled in ./external/psi3.tar.gz .

First, try
$ ./configure
$ make

If it works - great!

If the build fails, try to disable new features which rely on C++11 features and external packages
./configure  -DMPQC_NEW_FEATURES=OFF


Several configure options are available:
$ ./configure -h

Most relevant options/variables are:

--prefix - the installation directory
--debug - forces debug build
CXX - the C++ compiler
-D - passed directly to CMAKE command

The MPQC can either be built in source or in a build directory, e.g.
$ mdkir build
$ cd build
$ ../configure
$ make


For basic validation of MPQC do this:
$ make check0

For a thorougher validation of MPQC do this:
$ make check1

For the thoroughest validation of MPQC do this:
$ make check2

To install MPQC (program, libraries, and header files):
$ make install

To build and install programmer's documentation for MPQC:
$ cd doc
$ make install
