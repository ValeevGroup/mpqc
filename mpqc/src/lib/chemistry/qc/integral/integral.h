
#ifndef _chemistry_qc_integral_integral_h
#define _chemistry_qc_integral_integral_h

#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
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

  int has_side_effects();

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

  int has_side_effects();
  int has_side_effects_in_arg1();
  int has_side_effects_in_arg2();

  virtual ~OneBody3Int();
};

#endif
