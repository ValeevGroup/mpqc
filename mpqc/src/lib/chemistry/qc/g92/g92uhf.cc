
#include "g92.h"

////////////////////////////////////////////////////////////////////////////

#define CLASSNAME Gaussian92UHF
#define PARENTS public Gaussian92
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Gaussian92UHF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Gaussian92::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
Gaussian92UHF::init()
{
  estring_ = "E(UHF)";
  emethod_ = "#p units=au SCF=DIRECT UHF";
  gmethod_ = "#p units=au Force SCF=DIRECT UHF";
  hmethod_ = "#p units=au Freq SCF=DIRECT UHF";
}

Gaussian92UHF::Gaussian92UHF(const RefKeyVal&keyval):
  Gaussian92(keyval)
{
  if (!basis_) {
    fprintf(stderr,"Gaussian92UHF needs a basis\n");
    abort();
  }
  init();
}

Gaussian92UHF::~Gaussian92UHF()
{
}

Gaussian92UHF::Gaussian92UHF(StateIn&s) :
  Gaussian92(s)
  maybe_SavableState(s)
{
  init();
}

void
Gaussian92UHF::save_data_state(StateOut&s)
{
  Gaussian92::save_data_state(s);
}
