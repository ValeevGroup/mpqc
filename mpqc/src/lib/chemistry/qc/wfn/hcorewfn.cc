
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/integral/integralv2.h>

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

static void
occ(PointBag_double *z, int &nd, int &ns)
{
  int Z=0;
  for (Pix i=z->first(); i; z->next(i)) Z += (int) z->get(i);

  nd = Z/2;
  ns = Z%2;
}

HCoreWfn::HCoreWfn(StateIn& s) :
  OneBodyWavefunction(s)
{
  occ(_mol->charges(),ndocc,nsocc);
}

HCoreWfn::HCoreWfn(const RefKeyVal&keyval):
  OneBodyWavefunction(keyval)
{
  occ(_mol->charges(),ndocc,nsocc);
}

HCoreWfn::HCoreWfn(const OneBodyWavefunction& obwfn) :
  OneBodyWavefunction(obwfn)
{
  occ(_mol->charges(),ndocc,nsocc);
}

HCoreWfn::HCoreWfn(const HCoreWfn& hcwfn) :
  OneBodyWavefunction(hcwfn)
{
}

HCoreWfn::~HCoreWfn()
{
}

HCoreWfn &
HCoreWfn::operator=(const HCoreWfn& hcwfn)
{
  OneBodyWavefunction::operator=(hcwfn);
  ndocc = hcwfn.ndocc;
  nsocc = hcwfn.nsocc;
  return *this;
}

void
HCoreWfn::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
}

RefSCMatrix
HCoreWfn::eigenvectors()
{
  if (!_eigenvectors.computed()) {
    RefSymmSCMatrix h(basis_dimension());
    h.assign(0.0);
  
    // put T ints in h
    RefSCElementOp op = new GaussianKineticIntv2(basis(), molecule());
    h.element_op(op);
  
    // must do this since the new GaussianNuclearIntv2 is created before
    // the GaussianKineticIntv2 is destroyed, meaning that int_initialize_1e()
    // is called before int_done_1e().  This is a bad thing.
    op = 0;

    // calculate V ints
    RefSymmSCMatrix v(basis_dimension());
    v.assign(0.0);
  
    op = new GaussianNuclearIntv2(basis(), molecule());
    v.element_op(op);
    op=0;
  
    // Hcore = T+V
    h.accumulate(v);

    // free up memory
    v=0;
  
    // diagonalize Hcore, and transform to S^-1/2 basis
    RefSCMatrix vec(basis_dimension(), basis_dimension());
    RefDiagSCMatrix val(basis_dimension());
    h.diagonalize(val,vec);
  
    vec = ao_to_orthog_ao()*vec;

    _eigenvectors=vec;

    _eigenvectors.computed() = 1;
  }
  
  // this is a HACK!!!  For reasons I don't understand the GNU compiler is
  // not creating a new RefSCMatrix for the eigenvectors, so we have to
  // bump up the reference count ourselves.
  // _eigenvectors.result_noupdate()->reference();
  return _eigenvectors;
}

double
HCoreWfn::occupation(int i)
{

  if (i < ndocc)
    return 2.0;
  else if (i >= ndocc && i < ndocc+nsocc)
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
