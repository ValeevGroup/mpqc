
#ifndef _chemistry_qc_scf_ossscf_h
#define _chemistry_qc_scf_ossscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/scf.h>

////////////////////////////////////////////////////////////////////////////

class OSSSCF: public SCF {
#   define CLASSNAME OSSSCF
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
    int local_;

    ResultRefSymmSCMatrix cl_fock_;
    ResultRefSymmSCMatrix op_focka_;
    ResultRefSymmSCMatrix op_fockb_;

  public:
    OSSSCF(StateIn&);
    OSSSCF(const RefKeyVal&);
    ~OSSSCF();

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
    RefSymmSCMatrix cl_gmat_;
    RefSymmSCMatrix op_densa_;
    RefSymmSCMatrix op_densa_diff_;
    RefSymmSCMatrix op_gmata_;
    RefSymmSCMatrix op_densb_;
    RefSymmSCMatrix op_densb_diff_;
    RefSymmSCMatrix op_gmatb_;

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
SavableState_REF_dec(OSSSCF);

#endif
