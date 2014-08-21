# Description
This file describes how to compile MPQC.

# Want docs only?

To generate documentation without compiling the package do this in the source directory:  
$ bin/build_docs.sh  
Then open doc/html/index.html in a web browser.

# Prerequisites

MPQC is a complex piece of software and can be customized extensively. For example, to use MPQC on a parallel machine an MPI library
(see below) may be needed. While the basic subset of MPQC features will work on all platforms, advanced features of MPQC
will require additional prerequisites.

## Basic prerequisites

To compile even a bare-bones MPQC these prerequisites are needed.

* [CMake](http://www.cmake.org/)
  * CMake 2.8.8 or above is required. If no CMake is found in path, one will be downloaded and built in ./CMake when first running ./configure.
  * Running   
   $ ./configure cmake   
   will *force* building and using CMake in ./cmake directory
* Standard development tools; most of these should be available on any modern Unix-like system (Linux or OS X)
  * GNU Make
  * C/C++ and Fortran compilers
  * (optional) Message Passing Interface
  * (optional) POSIX Threads (Pthreads)
  * (optional) perl
  * (optional) python
* [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/)
  * It is preferable to use a vendor-provided optimized copy of BLAS+LAPACK,
    such as Intel Math Kernel Library (MKL) or Apple's Accelerate
    framework. A reference copy of BLAS+LAPACK will be downloaded and
    built only if not found on the system.

## Advanced prerequisites (optional)

More advanced features of MPQC that appeared for the first time in MPQC verion 3 require additional prerequisites. We recommend
that all of them are compiled with the same version of the C++ compiler as used to compile MPQC.

* C++ compiler that supports the C++11 standard; the most recent versions of gcc, clang, and icc are recommended.
* [Libint](https://github.com/evaleev/libint), version 2.1.0(beta) or later
  * Libint is bundled in ./external/libint*.tgz and will be built if requested. **NOTE** you may require custom version
  of Libint for some functionaility, such as interfacing with GAMESS.
* [TiledArray](https://github.com/ValeevGroup/tiledarray)
  * TiledArray must be compiled manually, with support for [Elemental](http://libelemental.org/)
* [Eigen](http://eigen.tuxfamily.org/)
  * A recent version is bundled in ./external/eigen.tar.gz; it will be unpacked if needed (this library is headers-only
  and does not need to be compiled).
* [Boost](http://www.boost.org/)
  * A recent version is bundled in ./external/boost*.tar.bz2; it will be built regardless of whether
    there is one is installed, since Boost depends heavily on compiler.
* [Psi3](http://www.psicode.org/)
  * Psi3 is bundled in ./external/psi3.tar.gz and will be compiled if needed. Psi4 does not work with MPQC.
* [HDF5](http://www.hdfgroup.org/HDF5/)
  * will be downloaded and built only if not found on the system.

# Build

First, try  
$ ./configure  
$ make  

If it works - great!

If the build fails, try to disable new features which rely on C++11 features and external packages:  
$ ./configure  -DMPQC_NEW_FEATURES=OFF

Several configure options are available:  
$ ./configure -h

Most relevant options/variables are:
* --prefix - the installation directory
* --debug - forces debug build
* CXX - the C++ compiler
* -D - passed directly to CMAKE command

The MPQC can either be built in source or in a build directory, e.g.  
$ mdkir build  
$ cd build  
$ ../configure  
$ make  

# Validate

For basic validation of MPQC do this:  
$ make check0

For a thorougher validation of MPQC do this:  
$ make check1

For the most thorough validation of MPQC do this:  
$ make check2

To install MPQC (program, libraries, and header files):  
$ make install

# Documentation (optional)

To build and install programmer's documentation for MPQC:  
$ cd doc  
$ make html install  
