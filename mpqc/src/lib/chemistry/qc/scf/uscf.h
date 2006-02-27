//
// uscf.h --- definition of the UnrestrictedSCF abstract base class
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_scf_uscf_h
#define _chemistry_qc_scf_uscf_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/scf.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/// A base class for unrestricted self-consistent-field methods.
class UnrestrictedSCF: public SCF {
  protected:
    Ref<PointGroup> most_recent_pg_;
    int user_occupations_;
    int tnalpha_;
    int tnbeta_;
    int nirrep_;
    int *nalpha_;
    int *nbeta_;
    int *initial_nalpha_;
    int *initial_nbeta_;

    AccResultRefSCMatrix oso_eigenvectors_beta_;
    AccResultRefDiagSCMatrix eigenvalues_beta_;
    ResultRefSymmSCMatrix focka_;
    ResultRefSymmSCMatrix fockb_;

  protected:
    Ref<SCExtrapError> extrap_error();
    // calculate the scf vector, returning the accuracy
    double compute_vector(double&, double enuclear);
    void initial_vector(int needv=1);
    
  public:
    UnrestrictedSCF(StateIn&);
    UnrestrictedSCF(const Ref<KeyVal>&);
    ~UnrestrictedSCF();

    void save_data_state(StateOut&);

    RefSCMatrix eigenvectors();
    RefDiagSCMatrix eigenvalues();

    RefSCMatrix oso_alpha_eigenvectors();
    RefSCMatrix alpha_eigenvectors();
    RefDiagSCMatrix alpha_eigenvalues();
    RefSCMatrix oso_beta_eigenvectors();
    RefSCMatrix beta_eigenvectors();
    RefDiagSCMatrix beta_eigenvalues();

    RefSymmSCMatrix alpha_density();
    RefSymmSCMatrix beta_density();
    RefSymmSCMatrix density();

    void symmetry_changed();

    double occupation(int, int);
    double alpha_occupation(int, int);
    double beta_occupation(int, int);
    
    // both return 1
    int spin_polarized();
    int spin_unrestricted();
    
    void print(std::ostream&o=ExEnv::out0()) const;

    int n_fock_matrices() const;
    /** Returns alpha (i==0) or beta (i==1) Fock matrix in AO basis (including XC contribution in KS DFT --
	compare this to CLSCF and HSOSSCF!). Argument i must be 0.
    */
    RefSymmSCMatrix fock(int i);
    /** Spin-unrestricted SCF methods do not define effective Fock matrix,
	thus this function should never be called. */
    RefSymmSCMatrix effective_fock();
    
    /** Overload of Function::set_desired_value_accuracy(). Must update
        accuracy of the eigenvalues and eigenvectors.
    */
    void set_desired_value_accuracy(double eps);
    
  protected:
    // these are temporary data, so they should not be checkpointed
    Ref<TwoBodyInt> tbi_;
    
    RefSymmSCMatrix densa_;
    RefSymmSCMatrix densb_;
    RefSymmSCMatrix gmata_;
    RefSymmSCMatrix gmatb_;
    RefSymmSCMatrix diff_densa_;
    RefSymmSCMatrix diff_densb_;

    void set_occupations(const RefDiagSCMatrix&);
    void set_occupations(const RefDiagSCMatrix&, const RefDiagSCMatrix&);

    void init_vector();
    void done_vector();
    double new_density();
    void reset_density();
    double scf_energy();
    Ref<SCExtrapData> extrap_data();

    void init_gradient();
    void done_gradient();
    RefSymmSCMatrix lagrangian();
    RefSymmSCMatrix gradient_density();
    
    void init_hessian();
    void done_hessian();

    // The Hartree-Fock derivatives
    void two_body_deriv_hf(double*grad,double exchange_fraction);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
