
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <util/keyval/keyval.h>
#include <math/nihmatrix/nihmatrix.h>
#include <math/nihmatrix/newmat.h>
#include <math/newmat7/cblas.h>
#include <math/newmat7/precisio.h>

#define CLASSNAME BFGSupdate
#define PARENTS public HessianUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BFGSupdate::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { HessianUpdate::_castdown(cd) };
  return do_castdowns(casts,cd);
}
BFGSupdate::BFGSupdate(KeyVal&kv):
  HessianUpdate(kv)
{
}
BFGSupdate::BFGSupdate(StateIn&s):
  SavableState(s,class_desc_),
  HessianUpdate(s)
{
  DVector x(s);
  Convert(x,xprev);
  DVector g(s);
  Convert(g,gprev);
}
void
BFGSupdate::save_data_state(StateOut&s)
{
  HessianUpdate::save_data_state(s);
  DVector tmp;
  Convert(xprev,tmp);
  tmp.save_object_state(s);
  Convert(gprev,tmp);
  tmp.save_object_state(s);
}

BFGSupdate::BFGSupdate()
{
}

BFGSupdate::~BFGSupdate()
{
}

BFGSupdate::BFGSupdate(SymmetricMatrix& H):
 HessianUpdate(H)
{
}

void
BFGSupdate::Initial(SymmetricMatrix& H,NLP1&nlp)
{
  xprev = nlp.GetXc();
  gprev = nlp.GetGrad();
  HessianUpdate::Initial(H,nlp);
}

void
BFGSupdate::Update(SymmetricMatrix& Hk,NLP1&nlp)
{
  Tracer trace("QNewton::UpdateH");

#ifdef DEBUG
  fprintf(stderr, "Inside BFGSupdate::UpdateH"\n);
  Print(sk);
#endif
  double mcheps = FloatingPointPrecision::Epsilon();

  int nr = nlp.GetDim();
  ColumnVector x(nlp.GetXc());
  ColumnVector grad(nlp.GetGrad());

// BFGS formula
  
  ColumnVector yk(nr), sk(nr), Bsk(nr);
  Matrix Htmp(nr,nr);
  
  yk = grad - gprev;
  sk = x   - xprev;
  
#ifdef DEBUG
  Print(yk);
  Print(sk);
#endif

  Real gts = Dot(gprev,sk);
  Real yts = Dot(yk,sk);
  
  Real snorm = Norm2(sk);
  Real ynorm = Norm2(yk);
  
#ifdef DEBUG
  printf("BFGSupdate: gts = %e, yts = %e\n",gts,yts);
  printf("BFGSupdate: snorm = %e, ynorm = %e\n",snorm,ynorm);
#endif

  if (yts <= sqrt(mcheps)*snorm*ynorm) {
    printf("BFGSupdate: <y,s> = %12.4e is too small\n", yts);
    printf("BFGSupdate: The BFGS update is skipped\n");
  }
  
  ColumnVector res(nr);
  res = yk - Hk*sk;
  if (res.NormInfinity() <= sqrt(mcheps)) {
    printf("BFGSupdate: <y,s> = %12.4e is too small\n", yts);
    printf("BFGSupdate: The BFGS update is skipped\n");
  }
  
// Otherwise update the Hessian approximation
#ifdef DEBUG
  printf("\nBFGSupdate: before update, k = %d\n", k);
  Print(Hk);
#endif

  Bsk = Hk*sk;
  Real sBs = Dot(sk,Bsk);
  
  Htmp = - (Bsk * Bsk.t()) / sBs;
  Htmp = Htmp + (yk * yk.t()) / yts;
  Htmp = Hk + Htmp;
  Hk << Htmp;
  Htmp.Release(); 
#ifdef DEBUG
  printf("\nBFGSupdate: after update, k = %d\n", k);
  Print(Hk);
  printf("BFGSupdate: sBs = %e\n",sBs);
#endif

  xprev = x;
  gprev = grad;
}
