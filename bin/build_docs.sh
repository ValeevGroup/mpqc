#!/bin/sh

function die(){
 echo $1
 exit 1
}

[ -e doc/Doxyfile.in ] || die "Did not find doc/Doxyfile.in, run this from the top-level MPQC source directory!"
command -v doxygen || die "Command \"doxygen\" is not found, make sure it's installed and in your PATH"

cd doc
sed s/@PROJECT_SOURCE_DIR@/../ < Doxyfile.in > Doxyfile.1
sed s/@MPQC_VERSION@/3.0.0-alpha/ < Doxyfile.1 > Doxyfile.2
rm -f Doxyfile.1
sed s/@HAVE_DOT@/NO/ < Doxyfile.2 > Doxyfile
rm -f Doxyfile.2
doxygen Doxyfile
rm -f Doxyfile
