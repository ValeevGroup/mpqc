
#ifndef _chemistry_qc_integral_integral_h
#define _chemistry_qc_integral_integral_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/integral/intiter.h>

class OneBodyInt: public SCElementOp
{
 protected:
  OneBodyIntIter *iter;
  RefGaussianBasisSet bs1;
  RefGaussianBasisSet bs2;
  double *buffer_;
 public:
  OneBodyInt(const RefGaussianBasisSet&b, OneBodyIntIter* =0);
  OneBodyInt(const RefGaussianBasisSet&b1, const RefGaussianBasisSet&b2,
             OneBodyIntIter* =0);
  
  virtual int nbasis();
  virtual int nbasis1();
  virtual int nbasis2();

  virtual int nshell();
  virtual int nshell1();
  virtual int nshell2();

  RefGaussianBasisSet basis() { return bs1; }
  RefGaussianBasisSet basis1() { return bs1; }
  RefGaussianBasisSet basis2() { return bs2; }

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
  OneBodyIntIter *iter;
  RefGaussianBasisSet bs1;
  RefGaussianBasisSet bs2;
  double *buffer_;
 public:
  OneBody3Int(const RefGaussianBasisSet&b, OneBodyIntIter* =0);
  OneBody3Int(const RefGaussianBasisSet&b1,const RefGaussianBasisSet&b2,
              OneBodyIntIter* =0);
  
  virtual int nbasis();
  virtual int nbasis1();
  virtual int nbasis2();
  virtual int nshell();
  virtual int nshell1();
  virtual int nshell2();

  RefGaussianBasisSet basis() { return bs1; }
  RefGaussianBasisSet basis1() { return bs1; }
  RefGaussianBasisSet basis2() { return bs2; }

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
