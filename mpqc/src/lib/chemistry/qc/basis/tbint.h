
#ifndef _chemistry_qc_basis_tbint_h
#define _chemistry_qc_basis_tbint_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <chemistry/qc/basis/gaussbas.h>

////////////////////////////////////////////////////////////////////////////

class TwoBodyInt : public VRefCount {
  protected:
    RefGaussianBasisSet bs1;
    RefGaussianBasisSet bs2;
    RefGaussianBasisSet bs3;
    RefGaussianBasisSet bs4;

    double *buffer_;

  public:
    TwoBodyInt(const RefGaussianBasisSet&b);
    TwoBodyInt(const RefGaussianBasisSet&b1,
               const RefGaussianBasisSet&b2,
               const RefGaussianBasisSet&b3,
               const RefGaussianBasisSet&b4);
    virtual ~TwoBodyInt();
  
    int nbasis() const;
    int nbasis1() const;
    int nbasis2() const;
    int nbasis3() const;
    int nbasis4() const;

    int nshell() const;
    int nshell1() const;
    int nshell2() const;
    int nshell3() const;
    int nshell4() const;

    RefGaussianBasisSet basis() { return bs1; }
    RefGaussianBasisSet basis1() { return bs1; }
    RefGaussianBasisSet basis2() { return bs2; }
    RefGaussianBasisSet basis3() { return bs3; }
    RefGaussianBasisSet basis4() { return bs4; }

    const double * buffer() const { return buffer_; }
    
    virtual void compute_shell(int,int,int,int) = 0;
};

REF_dec(TwoBodyInt);

////////////////////////////////////////////////////////////////////////////

class ShellQuartetIter {
  protected:
    const double * buf;
    double scale_;

    int e12;
    int e34;
    int e13e24;

    int index;
    
    int istart;
    int jstart;
    int kstart;
    int lstart;

    int iend;
    int jend;
    int kend;
    int lend;

    int icur;
    int jcur;
    int kcur;
    int lcur;

  public:
    ShellQuartetIter();
    virtual ~ShellQuartetIter();

    virtual void init(const double *,
                      int, int, int, int,
                      int, int, int, int,
                      int, int, int, int,
                      double);

    virtual void start();
    virtual void next();
    virtual operator int();

    int i() const;
    int j() const;
    int k() const;
    int l() const;

    int ij() const;
    int ik() const;
    int il() const;

    int kl() const;
    int jl() const;
    int jk() const;

    int ijkl() const;

    int nint() const { return iend*jend*kend*lend; }
    
    double val() const;
};

class TwoBodyIntIter {
  protected:
    RefTwoBodyInt tbi;
    ShellQuartetIter sqi;
    
    int iend;
    
    int icur;
    int jcur;
    int kcur;
    int lcur;
    
  public:
    TwoBodyIntIter();
    TwoBodyIntIter(const RefTwoBodyInt&);

    virtual ~TwoBodyIntIter();
    
    virtual void start();
    virtual void next();
    virtual operator int();

    int ishell() const;
    int jshell() const;
    int kshell() const;
    int lshell() const;

    virtual double scale() const;

    ShellQuartetIter& current_quartet();
};

#endif
