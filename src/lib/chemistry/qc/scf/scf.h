//
// scf.h --- definition of the SCF abstract base class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _chemistry_qc_scf_scf_h
#define _chemistry_qc_scf_scf_h

#include <util/group/thread.h>

#include <math/optimize/scextrap.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/wfn/accum.h>
#include <chemistry/qc/wfn/obwfn.h>

namespace sc {

  class DensityFittingInfo;

// //////////////////////////////////////////////////////////////////////////

/** The SCF class is the base for all classes that use a self-consistent
field procedure to solve an effective one body problem. */
class SCF: public OneBodyWavefunction {
  protected:
    int compute_guess_;

    int keep_guess_wfn_;
    Ref<OneBodyWavefunction> guess_wfn_;

    int always_use_guess_wfn_;

    Ref<SelfConsistentExtrapolation> extrap_;

    Ref<AccumH> accumdih_;
    Ref<AccumH> accumddh_;

    int maxiter_;
    int miniter_;
    int dens_reset_freq_;
    int reset_occ_;
    int local_dens_;
    size_t storage_;
    int print_all_evals_;
    int print_occ_evals_;

    double level_shift_;

    Ref<MessageGrp> scf_grp_;
    Ref<ThreadGrp> threadgrp_;
    int local_;

    Ref<TwoBodyInt>* tbis_; // a two body integral evaluator for each thread
    virtual void init_threads();
    virtual void done_threads();

    // implement the Compute::compute() function
    virtual void compute();

    // calculate the scf vector, returning the accuracy
    virtual double compute_vector(double&, double enuclear);

    // return the DIIS error matrices
    virtual Ref<SCExtrapError> extrap_error();

    // calculate the scf gradient
    virtual void compute_gradient(const RefSCVector&);

    // calculate the scf hessian
    virtual void compute_hessian(const RefSymmSCMatrix&);

    // saves state and restart information after every checkpoint_freq()
    // SCF iterations
    virtual void savestate_iter(int);

    // saves state to the given filename
    virtual void savestate_to_file(const std::string &filename);
    std::string previous_savestate_file_;

    // returns the log of the max density element in each shell block
    signed char * init_pmax(double *);

    // given a matrix, this will convert the matrix to a local matrix if
    // it isn't one already, and return that local matrix.  it will also
    // set the double* to point to the local matrix's data.
    enum Access { Read, Write, Accum };
    RefSymmSCMatrix get_local_data(const RefSymmSCMatrix&, double*&, Access);

    // create the initial scf vector.  either use the eigenvectors in
    // guess_wfn_, or use a core Hamiltonian guess.
    virtual void initial_vector();

    /// Obsolete scf vector so that next call to initial_vector() will cause recomputation
    virtual void obsolete_vector();

    // given the total number of density and fock matrices, figure out
    // how much memory that will require and then set the local_dens_
    // variable accordingly
    void init_mem(int);

    void so_density(const RefSymmSCMatrix& d, double occ, int alp=1);

    // Returns a new'ed allocation vector if it is in the input,
    // otherwise null.
    int *read_occ(const Ref<KeyVal> &, const char *name, int nirrep);

    /// how much lower is the desired accuracy of the guess?
    static double guess_acc_ratio() { return 1e4; }

    /// prints iteration log
    static void iter_print(int iter,
                           double energy,
                           double delta,
                           double walltime,
                           std::ostream& os = ExEnv::out0());

  public:
    SCF(StateIn&);
    /** The KeyVal constructor.

        <dl>

        <dt><tt>maxiter</tt><dd> This integer specifies the maximum number
        of SCF iterations.  The default is 40.

        <dt><tt>density_reset_frequency</tt><dd> This integer specifies how
        often, in term of SCF iterations, \f$\Delta D\f$ will be reset to
        \f$D\f$.  The default is 10.

        <dt><tt>reset_occupations</tt><dd> Reassign the occupations after
        each iteration based on the eigenvalues.  This only has an effect
        for molecules with higher than \f$C_1\f$ symmetry.  The default is
        false.

        <dt><tt>level_shift</tt><dd> Specificies a shift for the diagonal
        Fock matrix elements.  Doubly occupied orbitals are shifted by this
        amount and singly occupied orbitals are shifted by half this
        amount. This can improve convergence in difficult cases.  The units
        are hartrees and the default is 0.

        <dt><tt>extrap</tt><dd> This specifies an object of type
        SelfConsistentExtrapolation.  The default is a DIIS object.

        <dt><tt>memory</tt><dd> The amount of memory that each processor
        may use.  The default is 0 (minimal memory use).

        <dt><tt>local_density</tt><dd> If this is true, a local copy of the
        density and \f$G\f$ matrix will be made on all nodes, even if a
        distributed matrix specialization is used.  The default is true.

        <dt><tt>guess_wavefunction</tt><dd> This specifies the initial
        guess for the solution to the SCF equations.  This can be either a
        OneBodyWavefunction object or the name of file that contains the
        saved state of a OneBodyWavefunction object. By default,
        SuperpositionOfAtomicDensities object will be used; if that fails,
        core hamiltonian guess will be used instead.

        <dt><tt>keep_guess_wavefunction</tt><dd> The guess wavefunction is
        normally discarded after it is projected.  Setting this boolean
        variable to true will cause the guess to be kept.  This is useful
        when doing frequencies of symmetric molecules by finite
        displacements, because the wavefunction is lost whenever the
        molecule is displaced into lower symmetry.

        <dt><tt>always_use_guess_wavefunction</tt><dd> If the orbitals must
        be recomputed after they have already been computed once, then the
        old orbitals are used as the initial guess by default.  However, if
        this option is true, then the guess wavefunction will be used, if
        available.  If a guess wavefunction is not available, then a core
        Hamiltonian guess will be used.  If this option is set to true,
        then keep_guess_wavefunction should also be set to true.

        <dt><tt>print_evals</tt><dd>Takes a boolean value.  If true, print
        all eigenvalues after the SCF procedure converges.  Takes a boolean
        value.  The default is false.

        <dt><tt>print_occ_evals</tt><dd>Takes a boolean value.  If true,
        print the occupied eigenvalues after the SCF procedure converges.
        The default is false.

        <dt><tt>accumdih</tt><dd>Optional.  Takes an AccumH derivative.
        This provides additional contributions to the energy and the Fock
        matrix that are summed in once for the entire SCF procedure.
        AccumH's that are independent of the density should be given here.

        <dt><tt>accumddh</tt><dd>Optional.  Takes an AccumH derivative.
        This provides additional contributions to the energy and the Fock
        matrix that are summed in each iteration SCF procedure.  AccumH's
        that depend on the density must be given here.

        </dl> */
    SCF(const Ref<KeyVal>&);
    ~SCF();

    void save_data_state(StateOut&);

    RefSCMatrix oso_eigenvectors();
    RefDiagSCMatrix eigenvalues();

    int spin_unrestricted(); // return 0

    // return the number of AO Fock matrices needed
    virtual int n_fock_matrices() const =0;

    // returns the n'th AO Fock matrix
    virtual RefSymmSCMatrix fock(int) =0;

    // return the effective MO fock matrix
    virtual RefSymmSCMatrix effective_fock() =0;

    /// return the DensityFittingInfo object used to implement compute()
    /// this is important to be able to reconstruct the Fock matrix
    virtual Ref<DensityFittingInfo> dfinfo() const;

    virtual double one_body_energy();
    virtual void two_body_energy(double &ec, double &ex);

    void symmetry_changed();

    void obsolete();
    void purge();

    void print(std::ostream&o=ExEnv::out0()) const;

  protected:
    // the following are scratch and are not checkpointed
    RefSCMatrix oso_scf_vector_;
    RefSCMatrix oso_scf_vector_beta_; // only used if !spin_restricted
    RefSymmSCMatrix hcore_;

    // //////////////////////////////////////////////////////////////////////
    // pure virtual member functions follow

    // tries to automagically guess the MO occupations
    virtual void set_occupations(const RefDiagSCMatrix&) =0;

    // //////////////////////////////////////////////////////////////////////
    // do setup for SCF calculation
    virtual void init_vector() =0;
    virtual void done_vector() =0;

    // calculate new density matrices, returns the rms density difference
    virtual double new_density() =0;

    // reset density diff matrix and zero out delta G matrix
    virtual void reset_density() =0;

    // return the scf electronic energy
    virtual double scf_energy() =0;

    // return the initial extrapolation data. Used when the mixing_fraction
    // input for DIIS is nonzero. By default null is returned.
    virtual Ref<SCExtrapData> initial_extrap_data();

    // return the DIIS data matrices
    virtual Ref<SCExtrapData> extrap_data() =0;

    // form the AO basis fock matrices
    virtual void ao_fock(double accuracy) =0;

    // //////////////////////////////////////////////////////////////////////
    // do setup for gradient calculation
    virtual void init_gradient() =0;
    virtual void done_gradient() =0;

    virtual RefSymmSCMatrix lagrangian() =0;
    virtual RefSymmSCMatrix gradient_density() =0;
    virtual void two_body_deriv(double*) =0;

    // //////////////////////////////////////////////////////////////////////
    // do setup for hessian calculation
    virtual void init_hessian() =0;
    virtual void done_hessian() =0;

  private:
    // This experimental function does SVD of Coulomb matrix
    // to be used in low-rank reconstruction
    void svd_product_basis();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
