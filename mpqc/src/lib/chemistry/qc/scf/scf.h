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

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/optimize/scextrap.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/wfn/obwfn.h>

#define SCF_CHECK_BOUNDS 0

////////////////////////////////////////////////////////////////////////////

class SCF: public OneBodyWavefunction {
#   define CLASSNAME SCF
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefOneBodyWavefunction guess_wfn_;
    
    RefSCMatrix scf_vector_;
    RefSelfConsistentExtrapolation extrap_;
    
    int maxiter_;
    int dens_reset_freq_;
    int reset_occ_;
    int local_dens_;
    int storage_;
    
    double level_shift_;

    RefMessageGrp scf_grp_;
    int local_;

    int debug_;
    
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
    
    // create the initial scf vector.  either use the eigenvectors in
    // guess_wfn_, or use a core Hamiltonian guess.
    void initial_vector();
    
    // given the total number of density and fock matrices, figure out
    // how much memory that will require and then set the local_dens_
    // variable accordingly
    void init_mem(int);
    
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

    // nicely print n x 3 data that are stored in a vector
    void print_natom_3(const RefSCVector &, const char *t=0, ostream&o=cout);

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

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
