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

* Standard development tools; most of these should be available on any modern Unix-like system (Linux or OS X)
  * GNU Make
  * C/C++ and Fortran compilers
  * (optional) Message Passing Interface
  * (optional) POSIX Threads (Pthreads)
  * (optional) perl
  * (optional) python
* [CMake](http://www.cmake.org/)
  * CMake 2.8.8 or above is required. If no CMake is found in path, one will be downloaded and built in ./CMake when first running ./configure.
  * Running   
   $ ./configure cmake   
   will **force** building and using CMake in ./cmake directory
  there is one is installed, since Boost depends heavily on compiler.
* [LAPACK](http://www.netlib.org/lapack/)
  * will be downloaded and built only if not found on the system.

## Advanced prerequisites (optional)

More advanced features of MPQC that appeared for the first time in MPQC verion 3 require these additional packages:

* [TiledArray](https://github.com/ValeevGroup/tiledarray)
  * TiledArray must be compiled manually, with support for [Elemental](http://libelemental.org/)
* [Libint](https://github.com/evaleev/libint), version 2.0.5 or later
  * Libint is bundled in ./external/libint*.tgz and will be built if requested
* [Eigen](http://eigen.tuxfamily.org/)
  * A recent version is bundled in ./external/eigen.tar.gz; it will be unpacked if needed.
* [Boost](http://www.boost.org/)
  * A recent version is bundled in ./external/boost*.tar.bz2; it will be built regardless of whether
* [Psi3](http://www.psicode.org/)
  * Psi3 is bundled in ./external/psi3.tar.gz . Psi4 does not work with MPQC.
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
$ make install  
