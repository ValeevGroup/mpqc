
#include <iostream.h>

#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>

int
main(int argc, char *argv[])
{
  int i;

  char *infile = (argv[1]) ? argv[1] : "mpqc.in";
  RefKeyVal kv(new ParsedKeyVal(infile));

  RefMolecule mol = kv->describedclassvalue("molecule");

  mol->print();

  mol->move_to_com();
  cout << "Molecule at com:\n";
  mol->print();
  
  mol->transform_to_principal_axes(0);
  cout << "Molecule wrt principal axes:\n";
  mol->print();
  mol->point_group().symm_frame().print();

  mol->symmetrize();
  cout << "symmetrized molecule\n";
  mol->print();

  mol->cleanup_molecule();
  cout << "cleaned molecule\n";
  mol->print();
  
  int nunique = mol->num_unique_atoms();
  int * unique_atoms = mol->find_unique_atoms();

  cout << scprintf("\nnunique=%d: ",nunique);
  for (i=0; i < nunique; i++) cout << scprintf(" %d",unique_atoms[i]+1);
  cout << endl;

  RefMolecule unique = new Molecule;
  for (i=0; i < nunique; i++)
    unique->add_atom(i,mol->atom(unique_atoms[i]));

  cout << "unique atoms\n";
  unique->print();
  
  exit(0);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
