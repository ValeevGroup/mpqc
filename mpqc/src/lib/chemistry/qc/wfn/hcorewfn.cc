
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
  accumh(s)
  maybe_SavableState(s)
{
  occ(molecule()->charges(),ndocc,nsocc);
}

HCoreWfn::HCoreWfn(const RefKeyVal&keyval):
  OneBodyWavefunction(keyval)
{
  occ(molecule()->charges(),ndocc,nsocc);
  accumh = new AccumHCore;
  accumh->init(basis(),integral());
}

HCoreWfn::~HCoreWfn()
{
}

void
HCoreWfn::save_data_state(StateOut&s)
{
  OneBodyWavefunction::save_data_state(s);
  accumh.save_state(s);
}

RefSCMatrix
HCoreWfn::eigenvectors()
{
  if (!eigenvectors_.computed()) {

    // create the core Hamiltonian Hcore
    RefSymmSCMatrix h(basis_dimension(), basis_matrixkit());
    accumh->accum(h);
  
    // diagonalize Hcore, and transform to S^-1/2 basis
    RefSCMatrix vec(basis_dimension(), basis_dimension(), basis_matrixkit());
    RefDiagSCMatrix val(basis_dimension(), basis_matrixkit());
    h.diagonalize(val,vec);
  
    vec = ao_to_orthog_ao()*vec;

    eigenvectors_=vec;

    eigenvectors_.computed() = 1;
  }
  
  return eigenvectors_.result_noupdate();
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
