
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>
#include <util/keyval/keyval.h>
#include <math/newmat7/cblas.h>
#include <math/newmat7/precisio.h>

DescribedClass_REF_def(HessianUpdate);

#define CLASSNAME HessianUpdate
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
HessianUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}
HessianUpdate::HessianUpdate(KeyVal&kv)
{
  PrefixKeyVal pkv("guess",kv);
  if (pkv.exists(0,0)) {
      DMatrix tmp(pkv);
      guessH = new SymmetricMatrix;
      Convert(tmp,*guessH);
    }
  else {
      guessH = 0;
    }
}
void
HessianUpdate::save_data_state(StateOut&s)
{
  if (guessH) {
      s.put(1);
      DMatrix m;
      Convert(*guessH,m);
      m.save_object_state(s);
    }
  else {
      s.put(0);
    }
}
HessianUpdate::HessianUpdate(StateIn&s):
  SavableState(s,class_desc_)
{
  int have_guessH;
  s.get(have_guessH);
  if (have_guessH) {
      DMatrix h(s);
      guessH = new SymmetricMatrix;
      Convert(h,*guessH);
    }
  else {
      guessH = 0;
    }
}

HessianUpdate::HessianUpdate():
 guessH(0)
{
}

HessianUpdate::HessianUpdate(SymmetricMatrix&guess)
{
  guessH = new SymmetricMatrix(guess);
}

HessianUpdate::~HessianUpdate()
{
  if (guessH) delete guessH;
}

void HessianUpdate::Initial(SymmetricMatrix&H,NLP1&nlp)
{
  int nr = nlp.GetDim();
  if (guessH == 0) {
      // Initial Hessian is set equal to the Identity Matrix
      H.ReDimension(nr);
      H = 0.0;
      // The old hessian choice depended on sx which was just
      // a vector of 1.0's, so we still get the same H if sx
      // is omitted -- clj.
      //D = sx.AsDiagonal()*sx.AsDiagonal()*max(fabs(fvalue),1.0);
      for (int i=1; i <= nr; i++) H(i,i) = max(fabs(nlp.GetF()),1.0);
    }
  else {
      H = *guessH;
      delete guessH;
      guessH = 0;
    }
}
