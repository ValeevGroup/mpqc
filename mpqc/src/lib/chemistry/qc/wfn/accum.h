
#ifndef _chemistry_qc_wfn_accum_h
#define _chemistry_qc_wfn_accum_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/integral.h>

////////////////////////////////////////////////////////////////////////////

// computes the density independent part of H
class AccumDIH: public SavableState {
#   define CLASSNAME AccumDIH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefGaussianBasisSet basis_set_;
    RefIntegral integral_;

  public:
    AccumDIH();
    AccumDIH(StateIn&);
    AccumDIH(const RefKeyVal&);
    virtual ~AccumDIH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefGaussianBasisSet&, const RefIntegral&);
    virtual void accum(const RefSymmSCMatrix& h) =0;
    virtual void done();
};
SavableState_REF_dec(AccumDIH);

////////////////////////////////////////////////////////////////////////////

// computes the density dependent part of H
class AccumDDH: public SavableState {
#   define CLASSNAME AccumDDH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefGaussianBasisSet basis_set_;
    RefIntegral integral_;
    
  public:
    AccumDDH();
    AccumDDH(StateIn&);
    AccumDDH(const RefKeyVal&);
    virtual ~AccumDDH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefGaussianBasisSet&, const RefIntegral&);
    virtual void accum(const RefSymmSCMatrix& h,
                       const RefSymmSCMatrix& h_open) = 0;
    virtual void done();
};
SavableState_REF_dec(AccumDDH);

class AccumNullDDH: public AccumDDH {
#   define CLASSNAME AccumNullDDH
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    AccumNullDDH();
    AccumNullDDH(StateIn&);
    AccumNullDDH(const RefKeyVal&);
    ~AccumNullDDH();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h, const RefSymmSCMatrix& h_open);
};

SavableState_REF_dec(AccumNullDDH);

#endif
