

#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>
#include <chemistry/qc/scf/fock.h>

///////////////////////////////////////////////////////////////////////////
// Occupation

SavableState_REF_def(Occupation);

#define CLASSNAME Occupation
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#define PARENTS virtual_base public SavableState
#include <util/class/classi.h>
void *
Occupation::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Occupation::Occupation(const RefKeyVal& kv)
{
  n_=nd_=ns_=np_=0;
  nd=ns=np=0;
  occnum_=0;
  
  RefGaussianBasisSet gbs = kv->describedclassvalue("basis");
  init(gbs);
}

Occupation::Occupation(StateIn& s) :
  SavableState(s),
  occnum_(0)
{
  s.get(n_);
  s.get(nd_);
  s.get(ns_);
  s.get(np_);

  s.get(nd);
  s.get(ns);
  s.get(np);

  if (n_) {
    occnum_ = new double*[n_];
    for (int i=0; i < n_; i++)
      s.get(occnum_[i]);
  }
}

Occupation::~Occupation()
{
  if (nd) {
    delete[] nd;
    nd=0;
  }
  if (ns) {
    delete[] nd;
    nd=0;
  }
  if (np) {
    delete[] nd;
    nd=0;
  }

  if (occnum_) {
    for (int i=0; i < n_; i++)
      delete[] occnum_[i];
    delete[] occnum_;
    occnum_=0;
  }

  n_=nd_=ns_=np_=0;
}

void
Occupation::save_data_state(StateOut& s)
{
  s.put(n_);
  s.put(nd_);
  s.put(ns_);
  s.put(np_);

  if (n_) {
    s.put(nd,n_);
    s.put(ns,n_);
    s.put(np,n_);

    for (int i=0; i < n_; i++)
      s.put(occnum_[i],nd[i]+ns[i]+np[i]);
  }
}

void
Occupation::init(const RefGaussianBasisSet& gbs)
{
  SymmGaussianBasisSet *sgbs =
    SymmGaussianBasisSet::castdown(gbs);

  if (sgbs) {
  } else if (gbs.nonnull()) {
    n_=0;
  } else {
  }
}

double
Occupation::occupation(int i) const
{
  return 0;
}

double
Occupation::occupation(int ir, int i) const
{
  return 0;
}

int
Occupation::ndocc() const
{
  return 0;
}

int
Occupation::ndocc(int i) const
{
  return 0;
}

int
Occupation::nsocc() const
{
  return 0;
}

int
Occupation::nsocc(int i) const
{
  return 0;
}

int
Occupation::npocc() const
{
  return 0;
}

int
Occupation::npocc(int i) const
{
  return 0;
}
