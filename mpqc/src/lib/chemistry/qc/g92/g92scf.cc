
#include "g92.h"

////////////////////////////////////////////////////////////////////////////

static char *
new_string(const char *str)
{
  char *ret = new char[strlen(str)+1];
  strcpy(ret,str);
  return ret;
}

////////////////////////////////////////////////////////////////////////////

#define CLASSNAME Gaussian92SCF
#define PARENTS public Gaussian92
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Gaussian92SCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Gaussian92::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
Gaussian92SCF::init()
{
  estring_ = "E(RHF)";
  emethod_ = "#p units=au SCF=DIRECT RHF";
  gmethod_ = "#p units=au Force SCF=DIRECT RHF";
  hmethod_ = "#p units=au Freq SCF=DIRECT RHF";
}

Gaussian92SCF::Gaussian92SCF(const RefKeyVal&keyval):
  Gaussian92(keyval)
{
  if (!basis_) {
    fprintf(stderr,"Gaussian92SCF needs a basis\n");
    abort();
  }
  init();
}

Gaussian92SCF::~Gaussian92SCF()
{
}

Gaussian92SCF::Gaussian92SCF(StateIn&s) :
  Gaussian92(s)
  maybe_SavableState(s)
{
  init();
}

void
Gaussian92SCF::save_data_state(StateOut&s)
{
  Gaussian92::save_data_state(s);
}
