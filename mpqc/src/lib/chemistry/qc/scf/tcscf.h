
#ifndef _chemistry_qc_scf_tcscf_h
#define _chemistry_qc_scf_tcscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/scf.h>

////////////////////////////////////////////////////////////////////////////

class TCSCF: public SCF {
#   define CLASSNAME TCSCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
 protected:
    int user_occupations_;
    int tndocc_;
    int nirrep_;
    int *ndocc_;
    int osa_;
    int osb_;

    double occa_;
    double occb_;

    double ci1_;
    double ci2_;

    ResultRefSymmSCMatrix focka_;
    ResultRefSymmSCMatrix fockb_;
    ResultRefSymmSCMatrix ka_;
    ResultRefSymmSCMatrix kb_;
    
  public:
    TCSCF(StateIn&);
    TCSCF(const RefKeyVal&);
    ~TCSCF();

    void save_data_state(StateOut&);

    void print(ostream&o=cout);

    double occupation(int ir, int vectornum);

    int n_fock_matrices() const;
    RefSymmSCMatrix fock(int);
    RefSymmSCMatrix effective_fock();

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();

  protected:
    // these are temporary data, so they should not be checkpointed
    RefTwoBodyInt tbi_;

    RefSymmSCMatrix cl_dens_;
    RefSymmSCMatrix cl_dens_diff_;
    RefSymmSCMatrix op_densa_;
    RefSymmSCMatrix op_densa_diff_;
    RefSymmSCMatrix op_densb_;
    RefSymmSCMatrix op_densb_diff_;

    RefSymmSCMatrix ao_gmata_;
    RefSymmSCMatrix ao_gmatb_;
    RefSymmSCMatrix ao_ka_;
    RefSymmSCMatrix ao_kb_;
    
    RefSymmSCMatrix cl_hcore_;
    
    void set_occupations(const RefDiagSCMatrix& evals);

    // scf things
    void init_vector();
    void done_vector();
    void reset_density();
    double new_density();
    double scf_energy();

    void ao_fock();

    RefSCExtrapData extrap_data();
    
    // gradient things
    void init_gradient();
    void done_gradient();

    RefSymmSCMatrix lagrangian();
    RefSymmSCMatrix gradient_density();
    void two_body_deriv(double*);

    // hessian things

    void init_hessian();
    void done_hessian();
};
SavableState_REF_dec(TCSCF);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
