
#ifdef __GNUC__
#pragma implementation
#endif

#include <chemistry/qc/scf/scf.h>

///////////////////////////////////////////////////////////////////////////
// SCF

#define CLASSNAME SCF
#define PARENTS public OneBodyWavefunction
#include <util/class/classia.h>
void *
SCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCF::SCF(StateIn& s) :
  OneBodyWavefunction(s)
  maybe_SavableState(s)
{
  s.get(maxiter_);
}

SCF::SCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval)
{
  if (keyval->exists("maxiter"))
    maxiter_ = keyval->intvalue("maxiter");
}

SCF::~SCF()
{
}

void
SCF::save_data_state(StateOut& s)
{
  s.put(maxiter_);
}

RefSCMatrix
SCF::eigenvectors()
{
  return eigenvectors_.result();
}

void
SCF::print(ostream&o)
{
  OneBodyWavefunction::print(o);
}
