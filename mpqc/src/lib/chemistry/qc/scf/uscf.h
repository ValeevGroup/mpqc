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

////////////////////////////////////////////////////////////////////////////

class UnrestrictedSCF: public SCF {
#   define CLASSNAME UnrestrictedSCF
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    int user_occupations_;
    int tnalpha_;
    int tnbeta_;
    int nirrep_;
    int *nalpha_;
    int *nbeta_;
    
    AccResultRefSCMatrix cb_;
    AccResultRefDiagSCMatrix eb_;
    ResultRefSymmSCMatrix focka_;
    ResultRefSymmSCMatrix fockb_;

  protected:
    RefSCExtrapError extrap_error();
    void compute_vector(double&);
    void initial_vector(int needv=1);
    
  public:
    UnrestrictedSCF(StateIn&);
    UnrestrictedSCF(const RefKeyVal&);
    ~UnrestrictedSCF();

    void save_data_state(StateOut&);

    RefSCMatrix eigenvectors();
    RefDiagSCMatrix eigenvalues();

    RefSCMatrix alpha_eigenvectors();
    RefDiagSCMatrix alpha_eigenvalues();
    RefSCMatrix beta_eigenvectors();
    RefDiagSCMatrix beta_eigenvalues();

    RefSymmSCMatrix alpha_density();
    RefSymmSCMatrix beta_density();
    RefSymmSCMatrix density();
    
    double occupation(int, int);
    double alpha_occupation(int, int);
    double beta_occupation(int, int);
    
    // both return 1
    int spin_polarized();
    int spin_unrestricted();
    
    void print(ostream&o=cout);

    int n_fock_matrices() const;
    RefSymmSCMatrix fock(int);
    RefSymmSCMatrix effective_fock();
    
  protected:
    // these are temporary data, so they should not be checkpointed
    RefTwoBodyInt tbi_;
    
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
    RefSCExtrapData extrap_data();

    void init_gradient();
    void done_gradient();
    RefSymmSCMatrix lagrangian();
    RefSymmSCMatrix gradient_density();
    
    void init_hessian();
    void done_hessian();
};
SavableState_REF_dec(UnrestrictedSCF);

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
