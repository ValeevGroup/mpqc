
#include <math/scmat/local.h>
#include "molecule.h"
#include "energy.h"
#include "coor.h"

// force linkage of the taylor expansion energy evaluator
#include "taylor_f.h"

__builtin_delete(void*ptr)
{
  if (ptr>(void*)0 && ptr<(void*)0x100) abort();
  if (ptr) free(ptr);
}

void do_displacement(RefMolecularCoor&mc,int i);

main()
{
  RefKeyVal kv(new ParsedKeyVal(SRCDIR "/moltest.in"));
  RefKeyVal pkv(new PrefixKeyVal("molecule",kv));

  RefMolecule mol = new Molecule(pkv);

  mol_cleanup_molecule(*mol);
  printf("Clean Molecule:\n");
  mol->print();

  mol_transform_to_principal_axes(*mol);
  printf("Clean Molecule wrt principal axes:\n");
  mol->print();

  int nunique = mol_num_unique_atoms(*mol);
  int * unique_atoms = mol_find_unique_atoms(*mol);

  printf("nunique=%d: ",nunique);
  for (int i=0; i < nunique; i++) printf(" %d",unique_atoms[i]+1);
  printf("\n");

  mol->point_group().char_table().print();
exit(0);

  printf("getting simp:\n");
  RefKeyVal ppkv = new PrefixKeyVal("simp",kv);
  SetIntCoor simp(ppkv);
  simp.update_values(mol);
  printf("simp:\n");
  simp.print(mol);

  // compare the analytic bmatrix to the finite displacement bmatrix
  RefKeyVal pppkv = new PrefixKeyVal("bmat_test",kv);
  SetIntCoor bmat_test(pppkv);
  RefSCDimension dnc(new LocalSCDimension(bmat_test.n()));
  RefSCDimension dn3(new LocalSCDimension(mol->natom()*3));
  RefSCMatrix bmatrix(dnc,dn3);
  RefSCMatrix fd_bmatrix(dnc,dn3);
  printf("testing bmat with:\n");
  bmat_test.update_values(mol);
  bmat_test.print();
  bmat_test.bmat(mol,bmatrix);
  bmat_test.fd_bmat(mol,fd_bmatrix);
  cout << "test bmatrix:\n";
  bmatrix.print();
  cout << "difference between test and finite displacement bmatrix:\n";
  (fd_bmatrix - bmatrix).print();

  // now we get ambitious
  RefMolecularCoor mc = kv->describedclassvalue("molcoor");

  mc->print();

  do_displacement(mc,0);
  do_displacement(mc,1);
  do_displacement(mc,2);
  do_displacement(mc,3);

  RefSymmSCMatrix hessian(mc->dim());
  mc->guess_hessian(hessian);

  cout << "The guess hessian:\n";
  hessian.print();

  RefMolecularEnergy me = kv->describedclassvalue("energy");
  me->print();
}


void
do_displacement(RefMolecularCoor&mc,int i)
{
  // now try to displace the geometry
  RefSCVector internal(mc->dim());
  mc->to_internal(internal);
  cout << "The initial internal coordinates:\n";
  internal.print();
  internal(i) = internal(i) + 0.2;
  cout << "The new internal coordinates:\n";
  internal.print();
  mc->to_cartesian(internal);
  mc->to_internal(internal);
  cout << "The actual new internal coordinates:\n";
  internal.print();
}
