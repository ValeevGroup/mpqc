
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/wfn/accum.h>
#include <chemistry/qc/basis/integral.h>

///////////////////////////////////////////////////////////////////////////
// AccumDIH

#define CLASSNAME AccumDIH
#define PARENTS public SavableState
#include <util/class/classia.h>

void *
AccumDIH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumDIH::AccumDIH()
{
}

AccumDIH::AccumDIH(StateIn&s) :
  SavableState(s)
{
  basis_set_.restore_state(s);
  integral_.restore_state(s);
}

// for now I'm assuming that the specializations of AccumDIH will call
// init() so we won't read in the integral or basis set here
AccumDIH::AccumDIH(const RefKeyVal&)
{
}

AccumDIH::~AccumDIH()
{
}

void
AccumDIH::save_data_state(StateOut& s)
{
  basis_set_.save_state(s);
  integral_.save_state(s);
}

void
AccumDIH::init(const RefGaussianBasisSet& b, const RefIntegral& i)
{
  basis_set_ = b;
  integral_ = i;
}

void
AccumDIH::done()
{
}

///////////////////////////////////////////////////////////////////////////
// AccumDDH

#define CLASSNAME AccumDDH
#define PARENTS public SavableState
#include <util/class/classia.h>
void *
AccumDDH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumDDH::AccumDDH()
{
}

AccumDDH::AccumDDH(StateIn& s) :
  SavableState(s)
{
  basis_set_.restore_state(s);
  integral_.restore_state(s);
}

// for now I'm assuming that the specializations of AccumDDH will call
// init() so we won't read in the integral or basis set here
AccumDDH::AccumDDH(const RefKeyVal&)
{
}

AccumDDH::~AccumDDH()
{
}

void
AccumDDH::save_data_state(StateOut& s)
{
  basis_set_.save_state(s);
  integral_.save_state(s);
}

void
AccumDDH::init(const RefGaussianBasisSet& b, const RefIntegral& i)
{
  basis_set_ = b;
  integral_ = i;
}

void
AccumDDH::done()
{
}

///////////////////////////////////////////////////////////////////////////
// AccumDDH

#define CLASSNAME AccumNullDDH
#define PARENTS public AccumDDH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
AccumNullDDH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDDH::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumNullDDH::AccumNullDDH()
{
}

AccumNullDDH::AccumNullDDH(StateIn&s) :
  AccumDDH(s)
{
}

AccumNullDDH::AccumNullDDH(const RefKeyVal& keyval) :
  AccumDDH(keyval)
{
}

AccumNullDDH::~AccumNullDDH()
{
}

void
AccumNullDDH::save_data_state(StateOut& s)
{
  AccumDDH::save_data_state(s);
}

void
AccumNullDDH::accum(const RefSymmSCMatrix& h, const RefSymmSCMatrix& h_open)
{
}

///////////////////////////////////////////////////////////////////////////
// AccumEffectiveH

#define CLASSNAME AccumEffectiveH
#define PARENTS public SCElementOp2
#include <util/class/classia.h>
void *
AccumEffectiveH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp2::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumEffectiveH::AccumEffectiveH()
{
}

AccumEffectiveH::AccumEffectiveH(StateIn& s) :
  SCElementOp2(s)
{
  s.get(dbegin_);
  s.get(dfence_);
  s.get(sbegin_);
  s.get(sfence_);
}

AccumEffectiveH::AccumEffectiveH(const RefKeyVal& keyval)
{
#if 0
  char *name[2] = {"closed", "open"};

  for (int i=0; i<2; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<=j; k++) {
        coef_[index(i,j,k)] = keyval->doublevalue(name[i], j, k);
      }
    }
  }
#endif
}

AccumEffectiveH::~AccumEffectiveH()
{
}

int
AccumEffectiveH::index(int hindex, int shelli, int shellj)
{
  if (shellj > shelli) {
    int tmp = shelli;
    shelli = shellj;
    shellj = tmp;
  }

  return hindex * 9 + ((shelli+1)*shelli)/2 + shellj;
}

int
AccumEffectiveH::shell(int ibasis)
{
  if (dbegin_ <= ibasis && ibasis < dfence_) return 0;
  if (sbegin_ <= ibasis && ibasis < sfence_) return 1;
  return 2;
}

void
AccumEffectiveH::save_data_state(StateOut& s)
{
  SCElementOp2::save_data_state(s);
  s.put(dbegin_);
  s.put(dfence_);
  s.put(sbegin_);
  s.put(sfence_);
}

void
AccumEffectiveH::process(SCMatrixBlockIter&i,SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
    int ri = shell(i.i());
    int rj = shell(j.j());
    i.set(i.get() * coef_[index(0, ri, rj)]
          + j.get() * coef_[index(1, ri, rj)]);
  }
}

void
AccumEffectiveH::docc(int begin, int fence)
{
  dbegin_ = begin;
  dfence_ = fence;
}

void
AccumEffectiveH::socc(int begin, int fence)
{
  sbegin_ = begin;
  sfence_ = fence;
}
