
#include <stdio.h>
#include <chemistry/molecule/molecule.h>


int
main(int argc, char *argv[])
{
  char *infile = (argv[1]) ? argv[1] : "mpqc.in";
  RefKeyVal kv(new ParsedKeyVal(infile));

  RefMolecule mol = kv->describedclassvalue("molecule");

  printf("Molecule:\n");
  mol->print();

  mol_move_to_com(mol);
  printf("Molecule at com:\n");
  mol->print();
  
  mol_transform_to_principal_axes(mol,0);
  printf("Molecule wrt principal axes:\n");
  mol->print();
  mol->point_group().symm_frame().print();

  mol->symmetrize();
  printf("symmetrized molecule\n");
  mol->print();

  mol_cleanup_molecule(mol);
  printf("cleaned molecule\n");
  mol->print();
  
  int nunique = mol_num_unique_atoms(mol);
  int * unique_atoms = mol_find_unique_atoms(mol);

  printf("\nnunique=%d: ",nunique);
  for (int i=0; i < nunique; i++) printf(" %d",unique_atoms[i]+1);
  printf("\n");

  exit(0);
}
