
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
  _basis_set.restore_state(s);
  _molecule.restore_state(s);
}

// for now I'm assuming that the specializations of AccumDIH will call
// init() so we won't read in the molecule or basis set here
AccumDIH::AccumDIH(const RefKeyVal&)
{
}

AccumDIH::~AccumDIH()
{
}

void
AccumDIH::save_data_state(StateOut& s)
{
  _basis_set.save_state(s);
  _molecule.save_state(s);
}

void
AccumDIH::init(const RefGaussianBasisSet& b,
               const RefMolecule& m)
{
  _basis_set = b;
  _molecule = m;
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
  _basis_set.restore_state(s);
  _molecule.restore_state(s);
}

// for now I'm assuming that the specializations of AccumDDH will call
// init() so we won't read in the molecule or basis set here
AccumDDH::AccumDDH(const RefKeyVal&)
{
}

AccumDDH::~AccumDDH()
{
}

void
AccumDDH::save_data_state(StateOut& s)
{
  _basis_set.save_state(s);
  _molecule.save_state(s);
}

void
AccumDDH::init(const RefGaussianBasisSet& b,
               const RefMolecule& m)
{
  _basis_set = b;
  _molecule = m;
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
  s.get(_dbegin);
  s.get(_dfence);
  s.get(_sbegin);
  s.get(_sfence);
}

AccumEffectiveH::AccumEffectiveH(const RefKeyVal& keyval)
{
#if 0
  char *name[2] = {"closed", "open"};

  for (int i=0; i<2; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<=j; k++) {
        _coef[index(i,j,k)] = keyval->doublevalue(name[i], j, k);
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
  if (_dbegin <= ibasis && ibasis < _dfence) return 0;
  if (_sbegin <= ibasis && ibasis < _sfence) return 1;
  return 2;
}

void
AccumEffectiveH::save_data_state(StateOut& s)
{
  SCElementOp2::save_data_state(s);
  s.put(_dbegin);
  s.put(_dfence);
  s.put(_sbegin);
  s.put(_sfence);
}

void
AccumEffectiveH::process(SCMatrixBlockIter&i,SCMatrixBlockIter&j)
{
  for (i.reset(),j.reset(); i; ++i,++j) {
    int ri = shell(i.i());
    int rj = shell(j.j());
    i.set(i.get() * _coef[index(0, ri, rj)]
          + j.get() * _coef[index(1, ri, rj)]);
  }
}

void
AccumEffectiveH::docc(int begin, int fence)
{
  _dbegin = begin;
  _dfence = fence;
}

void
AccumEffectiveH::socc(int begin, int fence)
{
  _sbegin = begin;
  _sfence = fence;
}
