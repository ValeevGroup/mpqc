
#ifndef _chemistry_qc_scf_grscf_h
#define _chemistry_qc_scf_grscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/optimize/scextrap.h>
#include <chemistry/qc/wfn/obwfn.h>

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
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    double _coef[18];
    int _dbegin;
    int _dfence;
    int _sbegin;
    int _sfence;
    // hindex is 0 for the closed and 1 for the open shell fock matrix
    // shelli and shellj are 0 for closed, 1 for open, and 2 for virtual
    int index(int hindex, int shelli, int shellj);
    // converts a basis function number to a shell number
    int shell(int ibasis);
    double& coef(int i, int j, int k) { return _coef[index(i,j,k)]; }

  public:
    AccumEffectiveH(StateIn&);
    AccumEffectiveH(const RefKeyVal&);
    virtual ~AccumEffectiveH();

    void save_data_state(StateOut&);
    
    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&);

    void docc(int begin, int fence);
    void socc(int begin, int fence);
};
SavableState_REF_dec(AccumEffectiveH);

////////////////////////////////////////////////////////////////////////////

class GRSCF: public OneBodyWavefunction
{
#   define CLASSNAME GRSCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
    RefSCDimension  _orth_dim;
    RefSCDimension  _mo_dim;

    RefSelfConsistentExtrapolation _extrap;
    RefSCExtrapData _data;
    RefSCExtrapError _error;

    RefAccumDIH _accumdih;
    RefAccumDDH _accumddh;
    RefAccumEffectiveH _accumeffh;

    int _ndocc;
    int _nsocc;
    int _density_reset_freq;

    int _maxiter;
    int _eliminate;

    virtual void converge_eigenvectors();
    virtual double electronic_energy();
    virtual double nuclear_energy();
    virtual double rms_delta_density();
    virtual void form_density(int iter,
                              const RefSCMatrix& evec,
                              const RefSymmSCMatrix& P,
                              const RefSymmSCMatrix& P_open,
                              const RefSymmSCMatrix& DP,
                              const RefSymmSCMatrix& DP_open);
                              

    char *ckptdir;
    char *fname;
    FILE *outfile;

    void init();
    
  public:
    GRSCF(StateIn&);
    GRSCF(const RefKeyVal&);
    ~GRSCF();

    void save_data_state(StateOut&);

    void compute();
    
    RefSCMatrix eigenvectors();
    double occupation(int vectornum);

    void print(SCostream&o=SCostream::cout);
};
SavableState_REF_dec(GRSCF);

#endif
