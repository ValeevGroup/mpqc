
#ifndef _chemistry_qc_scf_tcscf_h
#define _chemistry_qc_scf_tcscf_h

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

class TCSCF: public OneBodyWavefunction
{
#   define CLASSNAME TCSCF
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

    RefSymmSCMatrix _focka;
    RefSymmSCMatrix _fockb;
    RefSymmSCMatrix _ka;
    RefSymmSCMatrix _kb;
    RefDiagSCMatrix _fock_evals;
    
    int _ndocc;

    double occa;
    double occb;

    double ci1;
    double ci2;

    int _density_reset_freq;

    int _maxiter;
    int _eliminate;

    char *ckptdir;
    char *fname;

    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix _gr_dens;
    RefSymmSCMatrix _gr_opa_dens;
    RefSymmSCMatrix _gr_opb_dens;
    RefSymmSCMatrix _gr_hcore;
    RefSCMatrix _gr_vector;
    
    void init();
    void compute();
    void do_vector(double&,double&);
    void form_ao_fock(centers_t *, double*, double&);
    void do_gradient(const RefSCVector&);
    
  public:
    TCSCF(StateIn&);
    TCSCF(const TCSCF&);
    TCSCF(const RefKeyVal&);
    TCSCF(const OneBodyWavefunction&);
    ~TCSCF();

    TCSCF& operator=(const TCSCF&);
    
    void save_data_state(StateOut&);

    void print(SCostream&o=SCostream::cout);

    RefSCMatrix eigenvectors();

    double occupation(int vectornum);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();
};
SavableState_REF_dec(TCSCF);

#endif
