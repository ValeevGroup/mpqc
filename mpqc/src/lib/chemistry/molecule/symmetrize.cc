
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

  mol_move_to_com(mol);
  cout << "Molecule at com:\n";
  mol->print();
  
  mol_transform_to_principal_axes(mol,0);
  cout << "Molecule wrt principal axes:\n";
  mol->print();
  mol->point_group().symm_frame().print();

  mol->symmetrize();
  cout << "symmetrized molecule\n";
  mol->print();

  mol_cleanup_molecule(mol);
  cout << "cleaned molecule\n";
  mol->print();
  
  int nunique = mol_num_unique_atoms(mol);
  int * unique_atoms = mol_find_unique_atoms(mol);

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
