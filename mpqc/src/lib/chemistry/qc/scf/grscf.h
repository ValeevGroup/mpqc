
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
#include <chemistry/qc/wfn/effh.h>

////////////////////////////////////////////////////////////////////////////

class GRSCF: public OneBodyWavefunction
{
#   define CLASSNAME GRSCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
    RefSelfConsistentExtrapolation _extrap;
    RefSCExtrapData _data;
    RefSCExtrapError _error;

    RefAccumDIH _accumdih;
    RefAccumDDH _accumddh;
    RefAccumEffectiveH _accumeffh;

    RefSymmSCMatrix _fock;
    RefSymmSCMatrix _op_fock;
    RefDiagSCMatrix _fock_evals;
    
    int _ndocc;
    int _nsocc;
    int _density_reset_freq;

    int _maxiter;
    int _eliminate;

    char *ckptdir;
    char *fname;

    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix _gr_dens;
    RefSymmSCMatrix _gr_dens_diff;
    RefSymmSCMatrix _gr_op_dens;
    RefSymmSCMatrix _gr_op_dens_diff;
    RefSymmSCMatrix _gr_gmat;
    RefSymmSCMatrix _gr_op_gmat;
    RefSymmSCMatrix _gr_hcore;
    RefSCMatrix _gr_vector;
    
    void init();
    virtual void compute();
    virtual void do_vector(double&,double&);
    virtual void form_ao_fock(double&);
    virtual double scf_energy();
    virtual void do_gradient(const RefSCVector&);
    
  public:
    GRSCF(StateIn&);
    GRSCF(const GRSCF&);
    GRSCF(const RefKeyVal&);
    GRSCF(const OneBodyWavefunction&);
    ~GRSCF();

    GRSCF& operator=(const GRSCF&);
    
    void save_data_state(StateOut&);

    void print(SCostream&o=SCostream::cout);

    RefSCMatrix eigenvectors();

    double occupation(int vectornum);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();
};
SavableState_REF_dec(GRSCF);

#endif
