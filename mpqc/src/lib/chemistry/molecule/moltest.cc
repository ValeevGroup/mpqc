
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/coor.h>
#include <util/render/object.h>
#include <util/render/oogl.h>

// force linkage of the taylor expansion energy evaluator
#include <chemistry/molecule/taylor_f.h>
// and of the renderer
#include <chemistry/molecule/molrender.h>

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

  RefMolecule mol = kv->describedclassvalue("molecule");

  mol_cleanup_molecule(mol);
  printf("Clean Molecule:\n");
  mol->print();

  mol_transform_to_principal_axes(mol);
  printf("Clean Molecule wrt principal axes:\n");
  mol->print();

  int nunique = mol_num_unique_atoms(mol);
  int * unique_atoms = mol_find_unique_atoms(mol);

  printf("nunique=%d: ",nunique);
  for (int i=0; i < nunique; i++) printf(" %d",unique_atoms[i]+1);
  printf("\n");

  mol->point_group().char_table().print();

  RefRender ren = kv->describedclassvalue("renderer");
  RefRenderedObject renmol = kv->describedclassvalue("renderedmolecule");
  if (ren.nonnull() && renmol.nonnull()) ren->render(renmol);

  //exit(0);

  printf("getting simp:\n");
  RefSetIntCoor simp = kv->describedclassvalue("simp");
  RefIntCoorGen gen = kv->describedclassvalue("generator");
  if (gen.nonnull()) {
      gen->print();
    }
  printf("simp before update:\n");
  simp->print(mol);
  simp->update_values(mol);
  printf("simp:\n");
  simp->print(mol);

  // compare the analytic bmatrix to the finite displacement bmatrix
  RefSetIntCoor bmat_test = kv->describedclassvalue("bmat_test");
  RefSCDimension dnc(new LocalSCDimension(bmat_test->n()));
  RefSCDimension dn3(new LocalSCDimension(mol->natom()*3));
  RefSCMatrix bmatrix(dnc,dn3);
  RefSCMatrix fd_bmatrix(dnc,dn3);
  printf("testing bmat with:\n");
  bmat_test->update_values(mol);
  bmat_test->print();
  bmat_test->bmat(mol,bmatrix);
  bmat_test->fd_bmat(mol,fd_bmatrix);
  cout << "test bmatrix:\n";
  bmatrix.print();
  cout << "fd bmatrix:\n";
  fd_bmatrix.print();
  cout << "difference between test and finite displacement bmatrix:\n";
  (fd_bmatrix - bmatrix).print();

  cout.flush();
  cerr.flush();
  fflush(stdout);
  fflush(stderr);
  
  // now we get ambitious
  RefMolecularCoor mc = kv->describedclassvalue("molcoor");
  cout.flush();
  cerr.flush();
  fflush(stdout);
  fflush(stderr);

  if (mc.nonnull()) {
      mc->print();

      cout.flush();
      cerr.flush();
      fflush(stdout);
      fflush(stderr);

      // do_displacement(mc,0);
      // do_displacement(mc,1);
      // do_displacement(mc,2);
      // do_displacement(mc,3);

      RefSymmSCMatrix hessian(mc->dim());
      mc->guess_hessian(hessian);

      // cout << "The guess hessian:\n";
      // hessian.print();
    }

  RefMolecularEnergy me = kv->describedclassvalue("energy");
  if (me.nonnull()) {
      me->print();
    }

  return 0;
}


void
do_displacement(RefMolecularCoor&mc,int i)
{
  if (i>=mc->dim().n()) return;
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
