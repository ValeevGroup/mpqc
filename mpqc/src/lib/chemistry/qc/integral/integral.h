
#ifndef _chemistry_qc_integral_integral_h
#define _chemistry_qc_integral_integral_h

#include <math/scmat/matrix.h>
#include <chemistry/qc/basis/basis.h>

class OneBodyInt: public SCElementOp
{
 private:
  RefGaussianBasisSet bs;
  double *buffer_;
 public:
  OneBodyInt(const RefGaussianBasisSet&b);
  
  virtual int nbasis();
  virtual int nshell();
  RefGaussianBasisSet basis() { return bs; }

  virtual void compute_shell(int,int,double*) = 0;

  virtual void process(SCMatrixBlockIter&);
  virtual void process(SCMatrixRectBlock*);
  virtual void process(SCMatrixLTriBlock*);

  virtual ~OneBodyInt();
};

class OneBody3Int: public SCElementOp3
{
 private:
  RefGaussianBasisSet bs;
  double *buffer_;
 public:
  OneBody3Int(const RefGaussianBasisSet&b);
  
  virtual int nbasis();
  virtual int nshell();
  RefGaussianBasisSet basis() { return bs; }

  virtual void compute_shell(int,int,double*) = 0;

  virtual void process(SCMatrixBlockIter&,
                       SCMatrixBlockIter&,
                       SCMatrixBlockIter&);
  virtual void process(SCMatrixRectBlock*,
                       SCMatrixRectBlock*,
                       SCMatrixRectBlock*);
  virtual void process(SCMatrixLTriBlock*,
                       SCMatrixLTriBlock*,
                       SCMatrixLTriBlock*);

  virtual ~OneBody3Int();
};

#endif
