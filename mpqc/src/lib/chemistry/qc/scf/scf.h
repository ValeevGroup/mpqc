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

#ifdef __GNUC__
#pragma interface
#endif

#include <util/group/thread.h>

#include <math/optimize/scextrap.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/wfn/accum.h>
#include <chemistry/qc/wfn/obwfn.h>

// //////////////////////////////////////////////////////////////////////////

/** The SCF class is the base for all classes that use a self-consistent
field procedure to solve an effective one body problem. */
class SCF: public OneBodyWavefunction {
#   define CLASSNAME SCF
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int need_vec_;
    int compute_guess_;

    RefOneBodyWavefunction guess_wfn_;
    
    RefSelfConsistentExtrapolation extrap_;
    
    RefAccumH accumdih_;
    RefAccumH accumddh_;
    
    int maxiter_;
    int dens_reset_freq_;
    int reset_occ_;
    int local_dens_;
    int storage_;
    int print_all_evals_;
    int print_occ_evals_;
    
    double level_shift_;

    RefMessageGrp scf_grp_;
    RefThreadGrp threadgrp_;
    int local_;
    
  protected:
    // implement the Compute::compute() function
    virtual void compute();

    // calculate the scf vector, returning the accuracy
    virtual double compute_vector(double&);

    // return the DIIS error matrices
    virtual RefSCExtrapError extrap_error();

    // calculate the scf gradient
    virtual void compute_gradient(const RefSCVector&);
    
    // calculate the scf hessian
    virtual void compute_hessian(const RefSymmSCMatrix&);
    
    // returns the log of the max density element in each shell block
    signed char * init_pmax(double *);
    
    // given a matrix, this will convert the matrix to a local matrix if
    // it isn't one already, and return that local matrix.  it will also
    // set the double* to point to the local matrix's data.
    enum Access { Read, Write, Accum };
    RefSymmSCMatrix get_local_data(const RefSymmSCMatrix&, double*&, Access);
    
    // create the initial scf vector.  either use the eigenvectors in
    // guess_wfn_, or use a core Hamiltonian guess.  Call this with needv
    // equal to 0 if you expect to call it twice with the same geometry
    // (eg. when calling from both set_occupations() and init_vector()).
    virtual void initial_vector(int needv=1);
    
    // given the total number of density and fock matrices, figure out
    // how much memory that will require and then set the local_dens_
    // variable accordingly
    void init_mem(int);
    
    void so_density(const RefSymmSCMatrix& d, double occ, int alp=1);

  public:
    SCF(StateIn&);
    /** @memo The KeyVal constructor.

        \begin{description}

        \item[maxiter] This integer specifies the maximum number of SCF
        iterations.  The default is 40.

        \item[density_reset_frequency] This integer specifies how often, in
        term of SCF iterations, $\Delta D$ will be reset to $D$.  The
        default is 10.

        \item[reset_occuptions] Reassign the occupations after each
        iteration based on the eigenvalues.  This only has an effect for
        molecules with higher than $C_1$ symmetry.  The default is false.

        \item[level_shift] The default is 0.

        \item[extrap] This specifies an object of type
        SelfConsistentExtrapolation.  The default is a DIIS object.

        \item[memory] The amount of memory that each processor may use.
        The default is 0 (minimal memory use).

        \item[local_density] If this is true, a local copy of the density
        and $G$ matrix will be made on all nodes, even if a distributed
        matrix specialization is used.  The default is true.

        \item[guess_wavefunction] This specifies the initial guess for the
        solution to the SCF equations.  This can be either a
        OneBodyWavefunction object or the name of file that contains the
        saved state of a OneBodyWavefunction object.  By default the
        one-electron hamiltonian will be diagonalized to obtain the initial
        guess.

        \end{description}
     */
    SCF(const RefKeyVal&);
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

    virtual double one_body_energy();
    virtual void two_body_energy(double &ec, double &ex);

    void print(ostream&o=cout) const;

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
    
    // return the DIIS data matrices
    virtual RefSCExtrapData extrap_data() =0;
    
    // form the AO basis fock matrices
    virtual void ao_fock() =0;

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
};
SavableState_REF_dec(SCF);

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
