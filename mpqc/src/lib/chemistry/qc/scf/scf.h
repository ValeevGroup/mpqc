
#ifndef _chemistry_qc_scf_scf_h
#define _chemistry_qc_scf_scf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/optimize/scextrap.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/wfn/obwfn.h>

////////////////////////////////////////////////////////////////////////////

class SCF: public OneBodyWavefunction {
#   define CLASSNAME SCF
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefSCMatrix scf_vector_;
    
    int maxiter_;
    int int_store_;
    int dens_reset_freq_;
    double level_shift_;

    RefMessageGrp scf_grp_;
    
  protected:
    // implement the Compute::compute() function
    virtual void compute();

    // calculate the scf vector
    virtual void compute_vector(double&);

    // return the DIIS error matrices
    virtual RefSCExtrapError extrap_error();

    // calculate the scf gradient
    virtual void compute_gradient(const RefSCVector&);
    
    // calculate the scf hessian
    virtual void compute_hessian(const RefSymmSCMatrix&);
    
    // returns the log of the max density element in each shell block
    char * init_pmax(double *);
    
    // given a matrix, this will convert the matrix to a local matrix if
    // it isn't one already, and return that local matrix.  it will also
    // set the double* to point to the local matrix's data.
    enum Access { Read, Write, Accum };
    RefSymmSCMatrix get_local_data(const RefSymmSCMatrix&, double*&, Access);
    
  public:
    SCF(StateIn&);
    SCF(const RefKeyVal&);
    ~SCF();

    void save_data_state(StateOut&);

    RefSCMatrix eigenvectors();

    // return the number of AO Fock matrices needed
    virtual int n_fock_matrices() const =0;

    // returns the n'th AO Fock matrix
    virtual RefSymmSCMatrix fock(int) =0;

    // return the effective MO fock matrix
    virtual RefSymmSCMatrix effective_fock() =0;
    
    void print(ostream&o=cout);

  protected:
    ////////////////////////////////////////////////////////////////////////
    // pure virtual member functions follow
    
    // tries to automagically guess the MO occupations
    virtual void set_occupations(const RefDiagSCMatrix&) =0;
    
    ////////////////////////////////////////////////////////////////////////
    // do setup for SCF calculation
    virtual void init_vector() =0;
    virtual void done_vector() =0;

    // calculate new density matrices, returns the rms density difference
    virtual double new_density() =0;

    // reset density diff matrix and zero out delta G matrix
    virtual void reset_density() =0;

    // return the scf electronic energy
    virtual double scf_energy() =0;
    
    // return the DIIS data matrices
    virtual RefSCExtrapData extrap_data() =0;
    
    // form the AO basis fock matrices
    virtual void ao_fock() =0;

    ////////////////////////////////////////////////////////////////////////
    // do setup for gradient calculation
    virtual void init_gradient() =0;
    virtual void done_gradient() =0;

    virtual RefSymmSCMatrix lagrangian() =0;
    virtual RefSymmSCMatrix gradient_density() =0;
    virtual void two_body_deriv(double*) =0;
    
    ////////////////////////////////////////////////////////////////////////
    // do setup for hessian calculation
    virtual void init_hessian() =0;
    virtual void done_hessian() =0;
};
SavableState_REF_dec(SCF);

#endif
