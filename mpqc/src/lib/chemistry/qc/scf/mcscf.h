
#ifndef _chemistry_qc_scf_mcscf_h
#define _chemistry_qc_scf_mcscf_h

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

class MCSCF: public OneBodyWavefunction
{
#   define CLASSNAME MCSCF
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

    RefSymmSCMatrix _fockc;
    RefSymmSCMatrix _focka;
    RefSymmSCMatrix _fockb;
    RefSymmSCMatrix _fockab;

    RefSymmSCMatrix _ka;
    RefSymmSCMatrix _kb;

    RefDiagSCMatrix _fock_evals;
    
    int _ndocc;
    int aorb;
    int borb;

    double occa;
    double occb;

    double ci1;
    double ci2;
    double ci3;

    int root;
    
    int _density_reset_freq;

    int _maxiter;
    int _eliminate;

    char *ckptdir;
    char *fname;

    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix _densc;
    RefSymmSCMatrix _densa;
    RefSymmSCMatrix _densb;
    RefSymmSCMatrix _densab;
    RefSymmSCMatrix _densab2;

    RefSymmSCMatrix _gr_hcore;
    RefSCMatrix _gr_vector;
    
    void init();
    void compute();
    void do_vector(double&,double&);
    void form_ao_fock(centers_t *, double*, double&);
    void do_gradient(const RefSCVector&);
    
  public:
    MCSCF(StateIn&);
    MCSCF(const MCSCF&);
    MCSCF(const RefKeyVal&);
    MCSCF(const OneBodyWavefunction&);
    ~MCSCF();

    MCSCF& operator=(const MCSCF&);
    
    void save_data_state(StateOut&);

    void print(SCostream&o=SCostream::cout);

    RefSCMatrix eigenvectors();

    double occupation(int vectornum);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();
};
SavableState_REF_dec(MCSCF);

#endif
