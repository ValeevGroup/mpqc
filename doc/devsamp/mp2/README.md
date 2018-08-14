This code example, written in C++, demonstrates how to augment the functionality of MPQC by adding a new class that implements MP2 energy.
This serves as an illustration of how to *plug extra functionality into MPQC*, by
extending an existing abstract class; in this case the MPQC LCAOWavefunction class is extended to
implement new quantum chemical methods. The new class can be used alongside
all other MPQC methods to, for example, compute forces by finite differences, etc.

1. To compile this example you must have compiled MPQC and installed it using `make install` command. This will install the code examples under `$prefix/doc/share/mpqc-${version}/examples directory`, where `$prefix` is the installation prefix provided to configure script during configuration of MPQC (the default is `/usr/local/mpqc/$mpqcversion`). The rest of instructions will assume that you are in the `mp2` subdirectory of that directory.
2. Configure the example by typing `cmake .` (for a clean restart remove file `CMakeCache.txt` and directory `CMakeFiles`); then compile the mp2 example by typing `make mp2`. This will create an executable file called `mp2`. To run it type `./mp2 ./mp2.json`.
3. To compare the MP2 energy against that obtained with MPQC, simply change object type from MP2 to RMP2 in `mp2.json`,
   change keyword "localize" to false (since RMP2 implements canonical MP2 energy only), and run the example again as `./mp2 ./mp2.json` (or you can run the main MPQC executable: `mpqc ./mp2.json`).
