
#include "molecule.h"
#include "energy.h"
#include <iostream.h>
#include <stdlib.h>

void try_energy(const char *);
void try_molecule(const char *);

int
main (int argc, char **argv)
{
	
  if (argc < 2) {
     cout << "usage: " << argv[0] << " <hostname>" << endl;
     exit (-1);
  }

  try_molecule(argv[1]);
  try_energy(argv[1]);
  return 0;
}

void
try_molecule(const char *host)
{
  C_Molecule* p;
  double natom;

  TRY {
    p = C_Molecule::_bind("testmol:Chemistry", host, IT_X);
  } CATCHANY {
    cerr << "Bind to object failed" << endl;
    cerr << "Unexpected exception " << IT_X << endl;
    exit(1);
  } ENDTRY

  TRY {
      p->keyval_create(
          "object<Molecule>: (\n"
          "  symmetry=c2v\n"
          "  {atoms geometry} = {\n"
          "    H  [   1.5  0.0   1.0 ]\n"
          "    O  [   0.0  0.0   0.0 ]\n"
          "    }\n"
          " )\n"
          , IT_X);
  } CATCHANY {
      cerr << "create of mol failed" << endl;
      cerr << "Unexpected exception " << IT_X << endl;
      exit(1);
  } ENDTRY

  TRY {
    // try to cout the atoms
    natom = p->natom(IT_X);
  } CATCHANY {
    // an error occurred while trying to read the energy
    cerr << "call to count the atoms failed" << endl;
    cerr << "Unexpected exception " << IT_X << endl;
    exit(1);
  } ENDTRY

  cout << "molecule had " << natom << " atoms" << endl;
}

void
try_energy(const char * host)
{
  C_MolecularEnergy* p;
  double e;
	
  TRY {
    p = C_MolecularEnergy::_bind("test:Chemistry", host, IT_X);
  } CATCHANY {
    cerr << "Bind to object failed" << endl;
    cerr << "Unexpected exception " << IT_X << endl;
    exit(1);
  } ENDTRY

  TRY {
      p->keyval_create(
          "object<MPSCF>: (\n"
          " molecule<Molecule>: (\n"
          "  symmetry=c2v\n"
          "  {atoms geometry} = {\n"
          "    H  [   1.5  0.0   1.0 ]\n"
          "    O  [   0.0  0.0   0.0 ]\n"
          "    }\n"
          "  )\n"
          " ckpt_freq = 100\n"
          " save_vector = no\n"
          " basis<GaussianBasisSet>: (\n"
          "    basisfiles = [ v2g90supp.ipv2 v2g90.ipv2 ]\n"
          "    molecule = $..:molecule\n"
          "    name = sto3g\n"
          "  )\n"
          " docc = 5\n"
          " local_P = yes\n"
          " restart = no\n"
          " )\n"
          , IT_X);
  } CATCHANY {
      cerr << "create of mole failed" << endl;
      cerr << "Unexpected exception " << IT_X << endl;
      exit(1);
  } ENDTRY

  TRY {
    // try to read the energy
    e = p->energy(IT_X);
  } CATCHANY {
    // an error occurred while trying to read the energy
    cerr << "call to read the energy failed" << endl;
    cerr << "Unexpected exception " << IT_X << endl;
    exit(1);
  } ENDTRY

  cout << "energy is " << e << endl;
}
