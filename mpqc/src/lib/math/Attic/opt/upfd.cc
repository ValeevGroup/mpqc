
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "opt.h"
#include <math/newmat7/cblas.h>
#include <math/newmat7/precisio.h>

FDupdate::FDupdate()
{
}

FDupdate::~FDupdate()
{
}

FDupdate::FDupdate(SymmetricMatrix& H):
 HessianUpdate(H)
{
}

void FDupdate::Initial(SymmetricMatrix& H,NLP1&nlp)
{
  Update(H,nlp);
}

void FDupdate::Update(SymmetricMatrix& Hk,NLP1& nlp)
{
  Tracer trace("FDNewton::UpdateH");

#ifdef DEBUG
  fprintf(stderr, "Inside OptFDNewton::UpdateH"\n);
  Print(Hk);
#endif
  Real mcheps = FloatingPointPrecision::Epsilon();
  Real sqrteps = sqrt(mcheps);

  int i;
  double hi;
  double xtmp;

  int nr = nlp.GetDim();

  ColumnVector gx(nr), grad(nr), xc;
  Matrix Htmp(nr,nr);
   
  xc = nlp.GetXc();
  gx = nlp.GetGrad();
  
  for (i=1; i<=nr; i++) {
    // scaling has been removed for now - clj
    //hi = sqrteps*max(fabs(xc(i)),sx(i));
    hi = sqrteps*fabs(xc(i));
    copysign(hi,xc(i));
    xtmp = xc(i);
    xc(i) = xtmp + hi;
    nlp.SetX(xc);
    grad = nlp.EvalG(); 
    Htmp.Column(i) << (grad - gx) / hi;
    xc(i) = xtmp;
  }
  
  nlp.SetX(xc);
  Hk << (Htmp.t() + Htmp)/2.0;

}
