
#include "g92.h"

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

char *
Gaussian92UHF::emethod()
{
  static char method[32];
  int conv = (int) -log10(desired_value_accuracy());

  sprintf(method,"uhf scf=direct scfcon=%d",conv);

  return method;
}

char *
Gaussian92UHF::gmethod()
{
  static char method[36];
  int conv = (int) -log10(desired_value_accuracy());
  
  sprintf(method,"uhf force scf=direct scfcon=%d",conv);

  return method;
}

char *
Gaussian92UHF::hmethod()
{
  static char method[48];
  int conv = (int) -log10(desired_value_accuracy());
  int hconv = (int) -log10(desired_hessian_accuracy());
  
  sprintf(method,"uhf freq scf=direct scfcon=%d cphf=conv=%d",conv,hconv);

  return method;
}

Gaussian92UHF::Gaussian92UHF(const RefKeyVal&keyval):
  Gaussian92(keyval)
{
  if (!basis_) {
    fprintf(stderr,"Gaussian92UHF needs a basis\n");
    abort();
  }
}

Gaussian92UHF::~Gaussian92UHF()
{
}

Gaussian92UHF::Gaussian92UHF(StateIn&s) :
  Gaussian92(s)
  maybe_SavableState(s)
{
}

void
Gaussian92UHF::save_data_state(StateOut&s)
{
  Gaussian92::save_data_state(s);
}
