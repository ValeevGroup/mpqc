
#include <stdio.h>
#include <util/misc/ieee.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molshape.h>
#include <chemistry/solvent/bem.h>

int
main()
{
  // Abort on floating point errors.
  ieee_trap_errors();

  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/bemtest.in");

  RefBEMSolvent solvent = keyval->describedclassvalue("solvent");

  solvent->init();
  solvent->init_system_matrix();
  solvent->done();

  solvent->init();
  solvent->init_system_matrix();
  solvent->done();

  ConnollyShape2::print_counts();
  CS2Sphere::print_counts();

  fflush(stdout);
  fflush(stderr);
  return 0;
}
