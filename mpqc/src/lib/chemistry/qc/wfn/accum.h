
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

////////////////////////////////////////////////////////////////////////////

// computes the density independent part of H
class AccumDIH: public SavableState {
#   define CLASSNAME AccumDIH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefGaussianBasisSet _basis_set;
    RefMolecule _molecule;
  public:
    AccumDIH();
    AccumDIH(StateIn&);
    AccumDIH(const RefKeyVal&);
    virtual ~AccumDIH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefGaussianBasisSet&, const RefMolecule&);
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
    RefGaussianBasisSet _basis_set;
    RefMolecule _molecule;
  public:
    AccumDDH();
    AccumDDH(StateIn&);
    AccumDDH(const RefKeyVal&);
    virtual ~AccumDDH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefGaussianBasisSet&, const RefMolecule&);
    virtual void accum(const RefSymmSCMatrix& h,
                       const RefSymmSCMatrix& h_open) = 0;
    virtual void done();
};
SavableState_REF_dec(AccumDDH);

////////////////////////////////////////////////////////////////////////////

class AccumEffectiveH: public SCElementOp2 {
#   define CLASSNAME AccumEffectiveH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    double _coef[18];
    int _dbegin;
    int _dfence;
    int _sbegin;
    int _sfence;

    virtual void init() =0;
    
    // hindex is 0 for the closed and 1 for the open shell fock matrix
    // shelli and shellj are 0 for closed, 1 for open, and 2 for virtual
    int index(int hindex, int shelli, int shellj);
    // converts a basis function number to a shell number
    int shell(int ibasis);

    double& coef(int i, int j, int k) { return _coef[index(i,j,k)]; }

  public:
    AccumEffectiveH();
    AccumEffectiveH(StateIn&);
    AccumEffectiveH(const RefKeyVal&);
    virtual ~AccumEffectiveH();

    void save_data_state(StateOut&);
    
    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&);

    void docc(int begin, int fence);
    void socc(int begin, int fence);
};
SavableState_REF_dec(AccumEffectiveH);


#endif
