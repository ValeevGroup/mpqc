
#ifndef _chemistry_qc_scf_fock_h
#define _chemistry_qc_scf_fock_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/matrix.h>
#include <math/optimize/scextrap.h>
#include <chemistry/qc/basis/symgaussbas.h>
#include <chemistry/qc/wfn/accum.h>

class Occupation : virtual_base public SavableState {
#   define CLASSNAME Occupation
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;

    int nd_;
    int ns_;
    int np_;

    int *nd;
    int *ns;
    int *np;
    
    double **occnum_;

  public:
    Occupation(const RefKeyVal&);
    Occupation(StateIn&);
    virtual ~Occupation();

    void save_data_state(StateOut&);

    void init(const RefGaussianBasisSet&);

    double occupation(int) const;
    double occupation(int,int) const;

    int ndocc() const;
    int ndocc(int) const;

    int nsocc() const;
    int nsocc(int) const;

    int npocc() const;
    int npocc(int) const;
};
SavableState_REF_dec(Occupation);

#if 0
class SCFDensity : virtual_base public SavableState {
#   define CLASSNAME SCFDensity
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    virtual void make_contribution() =0;
    
  public:
    SCFDensity(const RefKeyVal&);
    SCFDensity(StateIn&);
    virtual ~SCFDensity();
    
    void save_data_state(StateOut&);

    void reset(const RefSCMatrix& scfvec, const RefOccupation&);
};
SavableState_REF_dec(SCFDensity);

class Fock : public SavableState {
#   define CLASSNAME Fock
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefGaussianBasisSet gbs_;

    RefOccupation occ_;
    RefSCFDensity density_;
    
  public:
    Fock(const RefKeyVal&);
    Fock(StateIn&);
    virtual ~Fock();

    void save_data_state(StateOut&);

    int build(RefSCMatrix& scfvec, RefAccumDDH&);
    
    virtual RefSymmSCMatrix effective_fock() =0;

    virtual RefSCExtrapData extrap_data() =0;
    virtual RefSCExtrapError extrap_error() =0;
    
    virtual double energy(const RefSymmSCMatrix& dih) =0;
};
SavableState_REF_dec(Fock);
#endif
    
#endif
