
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/basis/integral.h>

/////////////////////////////////////////////////////////////////////////

static void
occ(PointBag_double *z, int &nd, int &ns)
{
  int Z=0;
  for (Pix i=z->first(); i; z->next(i)) Z += (int) z->get(i);

  nd = Z/2;
  ns = Z%2;
}

/////////////////////////////////////////////////////////////////////////

SavableState_REF_def(HCoreWfn);

#define CLASSNAME HCoreWfn
#define PARENTS public OneBodyWavefunction
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
HCoreWfn::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

HCoreWfn::HCoreWfn(StateIn& s) :
  OneBodyWavefunction(s),
  _accumh(s)
  maybe_SavableState(s)
{
  occ(_mol->charges(),ndocc,nsocc);
}

HCoreWfn::HCoreWfn(const RefKeyVal&keyval):
  OneBodyWavefunction(keyval)
{
  occ(_mol->charges(),ndocc,nsocc);
  _accumh = new AccumHCore;
  _accumh->init(basis(),molecule());
}

HCoreWfn::HCoreWfn(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  occ(_mol->charges(),ndocc,nsocc);
  _accumh = new AccumHCore;
  _accumh->init(basis(),molecule());
}

HCoreWfn::HCoreWfn(const HCoreWfn& hcwfn) :
  OneBodyWavefunction(hcwfn)
{
  occ(_mol->charges(),ndocc,nsocc);
  _accumh = hcwfn._accumh;
}

HCoreWfn::~HCoreWfn()
{
}

HCoreWfn &
HCoreWfn::operator=(const HCoreWfn& hcwfn)
{
  OneBodyWavefunction::operator=(hcwfn);
  _accumh = hcwfn._accumh;
  ndocc = hcwfn.ndocc;
  nsocc = hcwfn.nsocc;
  return *this;
}

void
HCoreWfn::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
  _accumh.save_state(s);
}

RefSCMatrix
HCoreWfn::eigenvectors()
{
  if (!_eigenvectors.computed()) {

    // create the core Hamiltonian Hcore
    RefSymmSCMatrix h(basis_dimension());
    _accumh->accum(h);
  
    // diagonalize Hcore, and transform to S^-1/2 basis
    RefSCMatrix vec(basis_dimension(), basis_dimension());
    RefDiagSCMatrix val(basis_dimension());
    h.diagonalize(val,vec);
  
    vec = ao_to_orthog_ao()*vec;

    _eigenvectors=vec;

    _eigenvectors.computed() = 1;
  }
  
  return _eigenvectors;
}

double
HCoreWfn::occupation(int i)
{

  if (i < ndocc)
    return 2.0;
  else if (i < ndocc+nsocc)
    return 1.0;
  else
    return 0.0;
}

void
HCoreWfn::compute()
{
  fprintf(stderr,"HCoreWfn::compute(): don't call me\n");
  abort();
}
