
#include "molei.h"
#include <iostream.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  MolE* p;
  double e;
	
  if (argc < 2) {
     cout << "usage: " << argv[0] << " <hostname>" << endl;
     exit (-1);
  }
  TRY {
    p = MolE::_bind("test", argv[1], IT_X);
  } CATCHANY {
    cerr << "Bind to object failed" << endl;
    cerr << "Unexpected exception " << IT_X << endl;
    exit(1);
  } ENDTRY

  TRY {
      p->create(
          "mole<MPSCF>: (\n"
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

  return 0;
}
