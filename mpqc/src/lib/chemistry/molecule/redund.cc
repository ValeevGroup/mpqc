
extern "C" {
#include <math.h>
};

#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/molecule/simple.h>

#include <math/topology/bitarray.h>

///////////////////////////////////////////////////////////////////////////
// members of RedundMolecularCoor

#define CLASSNAME RedundMolecularCoor
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public IntMolecularCoor
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
RedundMolecularCoor::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = IntMolecularCoor::_castdown(cd);
  return do_castdowns(casts,cd);
}

RedundMolecularCoor::RedundMolecularCoor(RefMolecule&mol):
  IntMolecularCoor(mol)
{
}

RedundMolecularCoor::RedundMolecularCoor(const RefKeyVal& keyval):
  IntMolecularCoor(keyval)
{
  init();
}

RedundMolecularCoor::RedundMolecularCoor(StateIn& s):
  IntMolecularCoor(s)
{
}

RedundMolecularCoor::~RedundMolecularCoor()
{
}

void
RedundMolecularCoor::save_data_state(StateOut&s)
{
  IntMolecularCoor::save_data_state(s);
}

void
RedundMolecularCoor::form_coordinates()
{
  int i;
  int nbonds = bonds_->n();
  int nbends = bends_->n();
  int ntors = tors_->n();
  int nouts = outs_->n();
  int nextras = extras_->n();
  for (i=0; i<nbonds; i++) all_->add(bonds_->coor(i));
  for (i=0; i<nbends; i++) all_->add(bends_->coor(i));
  for (i=0; i<ntors; i++) all_->add(tors_->coor(i));
  for (i=0; i<nouts; i++) all_->add(outs_->coor(i));
  for (i=0; i<nextras; i++) all_->add(extras_->coor(i));

  variable_ = all_;
  
  int nredundant = nbonds + nbends + ntors + nouts + nextras;
  int nfixed = fixed_->n();

  // see how many coords we expect
  int n3 = molecule_->natom()*3;
  int nunique = n3 - 6; // need to detect linear

  if (nredundant < nunique) {
    fprintf(stderr,"RedundMolecularCoor::form_coordinates: "
            "found too few redundant coordinates\n");
    fprintf(stderr,"  (the geometry is probably bad)\n");
    abort();
  }
  
  // put the fixed coordinates into the constant list
  nfixed = fixed_->n();
  for (i=0; i<nfixed; i++) {
    constant_->add(fixed_->coor(i));
  }
}

void
RedundMolecularCoor::guess_hessian(RefSymmSCMatrix&hessian)
{
  variable_->guess_hessian(molecule_,hessian);
}

RefSymmSCMatrix
RedundMolecularCoor::inverse_hessian(RefSymmSCMatrix& hessian)
{
  RefSCDimension dredun = hessian.dim();
  
  RefSCMatrix bmat(dredun,molecule_->dim_natom3());
  variable_->bmat(molecule_,bmat);
  
  RefSymmSCMatrix bmbt(dredun);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);
  RefSymmSCMatrix bmbt_i = bmbt.gi();

  bmat = dredun->create_matrix(dredun);
  bmat.assign(0.0);
  for (int i=0; i < dredun->n(); i++)
    bmat.set_element(i,i,1.0);
  bmat = bmbt * bmat;
  bmat = bmbt_i*bmat;
  return bmat * (bmat*hessian*bmat).gi() * bmat;
}
