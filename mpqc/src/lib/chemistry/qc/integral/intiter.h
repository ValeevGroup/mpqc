
#ifndef _chemistry_qc_integral_intiter_h
#define _chemistry_qc_integral_intiter_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/basis/symgaussbas.h>

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

class SymmOneBodyIntIter : public OneBodyIntIter {
  protected:
    PetiteList& pl;
    
  public:
    SymmOneBodyIntIter(const RefSymmGaussianBasisSet&);
    ~SymmOneBodyIntIter();

    void start();
    void next();

    void start_ltri();
    void next_ltri();

    double scale() const;
};

#endif
