
#ifndef _chemistry_qc_integral_integral_h
#define _chemistry_qc_integral_integral_h

#include <math/scmat/matrix.h>
#include <chemistry/qc/basis/basis.h>

class OneBodyInt: public SCSymmElementOp
{
 private:
  const GaussianBasisSet* bs;
  double *buffer_;
 public:
  OneBodyInt(const GaussianBasisSet*b);
  
  virtual int nbasis();
  virtual int nshell();
  virtual GaussianBasisSet& basis();

  virtual void compute_shell(int,int,double*) = 0;

  virtual void process(SCMatrixBlockIter&);
  virtual void process(SCMatrixRectBlock*);
  virtual void process(SCMatrixLTriBlock*);

  virtual ~OneBodyInt();
};

#endif
