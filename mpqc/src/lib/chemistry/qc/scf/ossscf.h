
#ifndef _chemistry_qc_scf_ossscf_h
#define _chemistry_qc_scf_ossscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/optimize/scextrap.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/wfn/effh.h>

#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

////////////////////////////////////////////////////////////////////////////

class OSSSCF: public OneBodyWavefunction
{
#   define CLASSNAME OSSSCF
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
    RefSymmSCMatrix _opa_fock;
    RefSymmSCMatrix _opb_fock;
    RefDiagSCMatrix _fock_evals;
    
    int _ndocc;
    int _density_reset_freq;

    int _maxiter;
    int _eliminate;

    char *ckptdir;
    char *fname;

    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix _gr_dens;
    RefSymmSCMatrix _gr_dens_diff;
    RefSymmSCMatrix _gr_gmat;
    RefSymmSCMatrix _gr_opa_dens;
    RefSymmSCMatrix _gr_opa_dens_diff;
    RefSymmSCMatrix _gr_opa_gmat;
    RefSymmSCMatrix _gr_opb_dens;
    RefSymmSCMatrix _gr_opb_dens_diff;
    RefSymmSCMatrix _gr_opb_gmat;
    RefSymmSCMatrix _gr_hcore;
    RefSCMatrix _gr_vector;
    
    void init();
    virtual void compute();
    virtual void do_vector(double&,double&);
    virtual void form_ao_fock(centers_t *, double*);
    virtual double scf_energy();
    virtual void do_gradient(const RefSCVector&);
    void form_density(const RefSCMatrix& vec,
                              const RefSymmSCMatrix& density,
                              const RefSymmSCMatrix& density_diff,
                              const RefSymmSCMatrix& open_densitya,
                              const RefSymmSCMatrix& open_densitya_diff,
                              const RefSymmSCMatrix& open_densityb,
                              const RefSymmSCMatrix& open_densityb_diff);
    
  public:
    OSSSCF(StateIn&);
    OSSSCF(const OSSSCF&);
    OSSSCF(const RefKeyVal&);
    OSSSCF(const OneBodyWavefunction&);
    ~OSSSCF();

    OSSSCF& operator=(const OSSSCF&);
    
    void save_data_state(StateOut&);

    void print(SCostream&o=SCostream::cout);

    RefSCMatrix eigenvectors();

    double occupation(int vectornum);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();
};
SavableState_REF_dec(OSSSCF);

#endif
