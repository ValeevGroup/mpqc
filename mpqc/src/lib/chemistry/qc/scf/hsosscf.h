
#ifndef _chemistry_qc_scf_hsosscf_h
#define _chemistry_qc_scf_hsosscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/scf.h>

////////////////////////////////////////////////////////////////////////////

class HSOSSCF: public SCF {
#   define CLASSNAME HSOSSCF
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int user_occupations_;
    int tndocc_;
    int tnsocc_;
    int nirrep_;
    int *ndocc_;
    int *nsocc_;
    int local_;

  public:
    HSOSSCF(StateIn&);
    HSOSSCF(const RefKeyVal&);
    ~HSOSSCF();

    void save_data_state(StateOut&);

    void print(ostream&o=cout);

    double occupation(int irrep, int vectornum);

    int value_implemented();
    int gradient_implemented();
    int hessian_implemented();

  protected:
    // these are temporary data, so they should not be checkpointed
    RefTwoBodyInt tbi_;
    
    RefSymmSCMatrix cl_fock_;
    RefSymmSCMatrix cl_dens_;
    RefSymmSCMatrix cl_dens_diff_;
    RefSymmSCMatrix cl_gmat_;
    RefSymmSCMatrix op_fock_;
    RefSymmSCMatrix op_dens_;
    RefSymmSCMatrix op_dens_diff_;
    RefSymmSCMatrix op_gmat_;

    RefSymmSCMatrix cl_hcore_;

    char *init_pmax(double*);
    
    void set_occupations(const RefDiagSCMatrix& evals);

    // scf things
    void init_vector();
    void done_vector();
    void reset_density();
    double new_density();
    double scf_energy();

    void ao_fock();

    RefSCExtrapData extrap_data();
    RefSymmSCMatrix effective_fock();
    
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
SavableState_REF_dec(HSOSSCF);

#endif
