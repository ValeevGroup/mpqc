
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
    int store1_;
    int store2_;
    int int_store_;
    
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

    RefGaussianBasisSet basis();
    RefGaussianBasisSet basis1();
    RefGaussianBasisSet basis2();
    RefGaussianBasisSet basis3();
    RefGaussianBasisSet basis4();

    const double * buffer() const;
    
    virtual void compute_shell(int,int,int,int) = 0;

    int int_store1(int i) { store1_=i; };
    int int_store2(int i) { store2_=i; };
    int integral_storage(int i) { int_store_=i; };
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

    int i_;
    int j_;
    int k_;
    int l_;
    
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

    int ready() const { return icur < iend; }

    int i() const { return i_; }
    int j() const { return j_; }
    int k() const { return k_; }
    int l() const { return l_; }

    int nint() const { return iend*jend*kend*lend; }
    
    double val() const { return buf[index]*scale_; }
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

    int ready() const { return (icur < iend); }

    int ishell() const { return icur; }
    int jshell() const { return jcur; }
    int kshell() const { return kcur; }
    int lshell() const { return lcur; }

    virtual double scale() const;

    ShellQuartetIter& current_quartet();
};

#endif
