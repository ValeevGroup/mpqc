
#ifndef _chemistry_qc_basis_obint_h
#define _chemistry_qc_basis_obint_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/elemop.h>

#include <chemistry/qc/basis/gaussbas.h>

////////////////////////////////////////////////////////////////////////////

class OneBodyInt : public VRefCount {
  protected:
    RefGaussianBasisSet bs1;
    RefGaussianBasisSet bs2;

    double *buffer_;

  public:
    OneBodyInt(const RefGaussianBasisSet&b);
    OneBodyInt(const RefGaussianBasisSet&b1,
               const RefGaussianBasisSet&b2);
    virtual ~OneBodyInt();
  
    int nbasis() const;
    int nbasis1() const;
    int nbasis2() const;

    int nshell() const;
    int nshell1() const;
    int nshell2() const;

    RefGaussianBasisSet basis() { return bs1; }
    RefGaussianBasisSet basis1() { return bs1; }
    RefGaussianBasisSet basis2() { return bs2; }

    const double * buffer() const { return buffer_; }
    
    virtual void compute_shell(int,int) = 0;
};

REF_dec(OneBodyInt);

////////////////////////////////////////////////////////////////////////////

class ShellPairIter {
  private:
    const double * buf;
    double scale_;

    int e12;

    int index;
    
    int ioffset;
    int joffset;

    int iend;
    int jend;

    int icur;
    int jcur;
    
  public:
    ShellPairIter();
    ~ShellPairIter();

    void init(const double * buffer, int ishell, int jshell,
              int ioff, int joff, int nfunci, int nfunci, int redund=0,
              double scale=1.0);

    void start() { icur=jcur=index=0; }
    int ready() const { return (icur < iend); }

    void next() {
      if (jcur < ((e12)?(icur):((jend)-1))) {
        index++;
        jcur++;
        return;
      }

      jcur=0;
      icur++;

      index = icur*jend;
    }

    int current_i() const { return icur; }
    int current_j() const { return jcur; }

    int i() const { return icur+ioffset; }
    int j() const { return jcur+joffset; }

    int nint() const { return iend*jend; }
    
    double val() const { return buf[index]*scale_; }
};

////////////////////////////////////////////////////////////////////////////

class OneBodyIntIter : public VRefCount {
  protected:
    RefOneBodyInt obi; // help me obi wan
    ShellPairIter spi;
    
    int redund;
    
    int istart;
    int jstart;
    
    int iend;
    int jend;

    int icur;
    int jcur;

    int ij;
    
  public:
    OneBodyIntIter();
    OneBodyIntIter(const RefOneBodyInt&);
    virtual ~OneBodyIntIter();
    
    virtual void start(int ist=0, int jst=0, int ien=0, int jen=0);
    virtual void next();

    int ready() const { return (icur < iend); }

    int ishell() const { return icur; }
    int jshell() const { return jcur; }

    int ijshell() const { return ij; }

    int redundant(int i) { int ret=redund; redund=i; return ret; }
    
    virtual double scale() const;

    RefOneBodyInt one_body_int() { return obi; }

    ShellPairIter& current_pair();
};

REF_dec(OneBodyIntIter);

////////////////////////////////////////////////////////////////////////////

class OneBodyIntOp: public SCElementOp {
  protected:
    RefOneBodyIntIter iter;

  public:
    OneBodyIntOp(const RefOneBodyInt&);
    OneBodyIntOp(const RefOneBodyIntIter&);
    virtual ~OneBodyIntOp();
  
    virtual void process(SCMatrixBlockIter&);
    virtual void process(SCMatrixRectBlock*);
    virtual void process(SCMatrixLTriBlock*);

    int has_side_effects();
};

class OneBody3IntOp: public SCElementOp3 {
  private:
    RefOneBodyIntIter iter;

  public:
    OneBody3IntOp(const RefOneBodyInt&b);
    OneBody3IntOp(const RefOneBodyIntIter&);
    virtual ~OneBody3IntOp();
  
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

};

#endif
