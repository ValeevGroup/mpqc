
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/wfn/hcore.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>

#define CLASSNAME AccumHCore
#define PARENTS public AccumDIH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
AccumHCore::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDIH::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumHCore::AccumHCore()
{
}

AccumHCore::AccumHCore(StateIn&s) :
  AccumDIH(s)
{
}

AccumHCore::AccumHCore(const RefKeyVal& keyval) :
  AccumDIH(keyval)
{
}

AccumHCore::~AccumHCore()
{
}

void
AccumHCore::save_data_state(StateOut& s)
{
  AccumDIH::save_data_state(s);
}

void
AccumHCore::accum(const RefSymmSCMatrix& h)
{
  integral_->set_basis(basis_set_);

  RefSCElementOp hc = new OneBodyIntOp(integral_->kinetic());
  h.assign(0.0);
  h.element_op(hc);
  hc=0;

  RefOneBodyInt nuc = integral_->nuclear();
  nuc->reinitialize();
  hc = new OneBodyIntOp(nuc);
  h.element_op(hc);
  hc=0;
}

//////////////////////////////////////////////////////////////////////////////

#define CLASSNAME SymmAccumHCore
#define PARENTS public AccumDIH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
SymmAccumHCore::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumDIH::_castdown(cd);
  return do_castdowns(casts,cd);
}

SymmAccumHCore::SymmAccumHCore()
{
}

SymmAccumHCore::SymmAccumHCore(StateIn&s) :
  AccumDIH(s)
{
}

SymmAccumHCore::SymmAccumHCore(const RefKeyVal& keyval) :
  AccumDIH(keyval)
{
}

SymmAccumHCore::~SymmAccumHCore()
{
}

void
SymmAccumHCore::save_data_state(StateOut& s)
{
  AccumDIH::save_data_state(s);
}

void
SymmAccumHCore::accum(const RefSymmSCMatrix& h)
{
  integral_->set_basis(basis_set_);
  RefPetiteList pl = integral_->petite_list();

  // form skeleton Hcore in AO basis
  RefSymmSCMatrix hao(basis_set_->basisdim(), basis_set_->matrixkit());
  hao.assign(0.0);

  RefSCElementOp hc =
    new OneBodyIntOp(new SymmOneBodyIntIter(integral_->kinetic(), pl));
  hao.element_op(hc);
  hc=0;

  RefOneBodyInt nuc = integral_->nuclear();
  nuc->reinitialize();
  hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
  hao.element_op(hc);
  hc=0;

  // now symmetrize Hao
  pl->symmetrize(hao,h);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
