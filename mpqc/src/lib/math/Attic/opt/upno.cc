
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <util/keyval/keyval.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>
#include <math/newmat7/cblas.h>
#include <math/newmat7/precisio.h>

#define CLASSNAME NoUpdate
#define PARENTS public HessianUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
NoUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { HessianUpdate::_castdown(cd) };
  return do_castdowns(casts,cd);
}
NoUpdate::NoUpdate(KeyVal&kv):
  HessianUpdate(kv)
{
}
NoUpdate::NoUpdate(StateIn&s):
  SavableState(s,class_desc_),
  HessianUpdate(s)
{
}
void
NoUpdate::save_data_state(StateOut&s)
{
  HessianUpdate::save_data_state(s);
}

NoUpdate::NoUpdate()
{
}

NoUpdate::~NoUpdate()
{
}

NoUpdate::NoUpdate(SymmetricMatrix& H):
 HessianUpdate(H)
{
}

void
NoUpdate::Update(SymmetricMatrix& Hk,NLP1&nlp)
{
}
