
#ifndef _chemistry_qc_basis_obint_h
#define _chemistry_qc_basis_obint_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>

#include <chemistry/qc/basis/gaussbas.h>

////////////////////////////////////////////////////////////////////////////

class OneBodyIntIter {
  protected:
    int istart;
    int iend;
    int jstart;
    int jend;

    int icur;
    int jcur;
    
  public:
    OneBodyIntIter();
    OneBodyIntIter(int,int,int,int);
    virtual ~OneBodyIntIter();
    
    virtual void reset(int,int,int,int);
    
    virtual void start();
    virtual int ready();
    virtual void next();

    virtual void start_ltri();
    virtual int ready_ltri();
    virtual void next_ltri();

    virtual int ishell() const;
    virtual int jshell() const;

    virtual double scale() const;
};

////////////////////////////////////////////////////////////////////////////

class OneBodyInt: public SCElementOp {
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

class OneBody3Int: public SCElementOp3 {
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
