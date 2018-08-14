This code example, written in C++, demonstrates how to implement a simple standalone program, with its own `main()` function, by using the MPQC infrastructure.
This is an illustration of how to *plug MPQC into an existing codebase*;
specificaly, MPQC AO Integral Factories are used to compute and manipulate AO integrals.

1. To compile this example you must have compiled MPQC and installed it using `make install` command. This will install the code examples under `$prefix/doc/share/mpqc-${version}/examples directory`, where `$prefix` is the installation prefix provided to configure script during configuration of MPQC (the default is `/usr/local/mpqc/$mpqcversion`). The rest of instructions will assume that you are in the `aoints` subdirectory of that directory.
2. Configure the example by typing `cmake .`; then compile the aoints example by typing `make`. This will create an executable file called `aoints`. To run it type `./aoints ./aoints.json`.
