
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
  init();
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
  variable_ = all_;
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
  
  // form bmat for variable coordinates (ie all the simples)
  RefSCMatrix bmat(dredun,dnatom3_);
  variable_->bmat(molecule_,bmat);
  
  // and form G = (B*B+)
  RefSymmSCMatrix bmbt(dredun);
  bmbt.assign(0.0);
  bmbt.accumulate_symmetric_product(bmat);

  // free bmat, and allocate storage for the projection matrix p
  bmat = 0;

  RefSCMatrix p(dredun,dredun);
  p.assign(0.0);

  // form p = G- * G
  for (int i=0; i < dredun->n(); i++)
    p.set_element(i,i,1.0);
  p = bmbt * p;
  p = bmbt.gi()*p;
  
  // accumulate (p*hessian*p).gi() into bmbt
  bmbt.assign(0.0);
  bmbt.accumulate_transform(p,hessian);
  bmbt = bmbt.gi();
  
  // finally return hinv = p*(p*h*p)-*p
  RefSymmSCMatrix thess = hessian.clone();
  thess.assign(0.0);
  thess.accumulate_transform(p,bmbt);
  return thess;
}
