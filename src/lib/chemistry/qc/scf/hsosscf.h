//
// hsosscf.h --- definition of the high-spin open shell SCF class
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

#ifndef _chemistry_qc_scf_hsosscf_h
#define _chemistry_qc_scf_hsosscf_h

#include <chemistry/qc/scf/scf.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/** The HSOSSCF class is a base for classes implementing a self-consistent
procedure for high-spin open-shell molecules. */
class HSOSSCF: public SCF {
  protected:
    Ref<PointGroup> most_recent_pg_;
    int user_occupations_;
    int tndocc_;
    int tnsocc_;
    int nirrep_;
    int *initial_ndocc_;
    int *initial_nsocc_;
    int *ndocc_;
    int *nsocc_;

    ResultRefSymmSCMatrix cl_fock_;
    ResultRefSymmSCMatrix op_fock_;

    // computes semicanonical orbitals. Uses fock() to get the fock matrix.
    // Should be overloaded if fock() does not return the correct matrix (e.g. HSOSKS)
    virtual void semicanonical();
    bool semicanonical_evaluated_;
    AccResultRefDiagSCMatrix alpha_semican_evals_;
    AccResultRefDiagSCMatrix beta_semican_evals_;
    AccResultRefSCMatrix alpha_semican_evecs_;
    AccResultRefSCMatrix beta_semican_evecs_;

  public:
    HSOSSCF(StateIn&);
    /** The KeyVal constructor.
        <dl>

        <dt><tt>total_charge</tt><dd> This floating point number
        gives the total charge,
        \f$c\f$, of the molecule.  The default is 0.

        <dt><tt>nsocc</tt><dd> This integer gives the total number of
        singly occupied orbitals, \f$n_\mathrm{socc}\f$.  If this is not
        given, then multiplicity will be read.

        <dt><tt>multiplicity</tt><dd> This integer gives the multiplicity,
        \f$m\f$, of the molecule.  The number of singly occupied orbitals
        is then \f$n_\mathrm{socc} = m - 1\f$.  If neither nsocc nor
        multiplicity is given, then the maximum value of multiplicity
        will be sought that minimizes the total sum of orbital energies from guess_wavefunction
        (multiplicity = 1 will not be considered in the search). \sa HundsFEMOSeeker

        <dt><tt>ndocc</tt><dd> This integer gives the total number of
        doubly occupied orbitals \f$n_\mathrm{docc}\f$.  The default
        \f$n_\mathrm{docc} = (c - n_\mathrm{socc})/2\f$.

        <dt><tt>socc</tt><dd> This vector of integers gives the total
        number of singly occupied orbitals of each irreducible
        representation.  By default, the \f$n_\mathrm{socc}\f$ singly
        occupied orbitals will be distributed according to orbital
        eigenvalues.  If socc is given, then docc must be given and they
        override nsocc, multiplicity, ndocc, and total_charge.

        <dt><tt>docc</tt><dd> This vector of integers gives the total
        number of doubly occupied orbitals of each irreducible
        representation.  By default, the \f$n_\mathrm{docc}\f$ singly
        occupied orbitals will be distributed according to orbital
        eigenvalues.  If docc is given, then socc must be given and they
        override nsocc, multiplicity, ndocc, and total_charge.

        <dt><tt>maxiter</tt><dd> This has the same meaning as in the parent
        class, SCF; however, the default value is 100.

        <dt><tt>level_shift</tt><dd> This has the same meaning as in the
        parent class, SCF; however, the default value is 1.0.

        </dl> */
    HSOSSCF(const Ref<KeyVal>&);
    ~HSOSSCF();

    void save_data_state(StateOut&);

    void print(std::ostream&o=ExEnv::out0()) const;

    double occupation(int irrep, int vectornum);
    double alpha_occupation(int irrep, int vectornum);
    double beta_occupation(int irrep, int vectornum);

    int n_fock_matrices() const;
    /** Returns closed-shell (i==0) or open-shell (i==1) Fock matrix in SO basis
        (excluding XC contribution in KS DFT).
	Use effective_fock() if you want the full KS Fock matrix.
    */
    RefSymmSCMatrix fock(int i);
    /** Returns effective Fock matrix in MO basis (including XC contribution for KS DFT). */
    RefSymmSCMatrix effective_fock();

    void symmetry_changed();
    
    // returns 1
    int spin_polarized();
    RefSymmSCMatrix density();
    RefSymmSCMatrix alpha_density();
    RefSymmSCMatrix beta_density();
    
    /// Coefficients of semicanonical alpha-spin orbitals in SO basis
    RefSCMatrix alpha_semicanonical_eigenvectors();
    /// Coefficients of semicanonical beta-spin orbitals in SO basis
    RefSCMatrix beta_semicanonical_eigenvectors();
    /// Eigenvalues of semicanonical alpha-spin orbitals
    RefDiagSCMatrix alpha_semicanonical_eigenvalues();
    /// Eigenvalues of semicanonical beta-spin orbitals
    RefDiagSCMatrix beta_semicanonical_eigenvalues();

  protected:
    // these are temporary data, so they should not be checkpointed
    RefSymmSCMatrix cl_dens_;
    RefSymmSCMatrix cl_dens_diff_;
    RefSymmSCMatrix cl_gmat_;
    RefSymmSCMatrix op_dens_;
    RefSymmSCMatrix op_dens_diff_;
    RefSymmSCMatrix op_gmat_;

    RefSymmSCMatrix cl_hcore_;

    /** Implementation of SCF::set_occupations(const RefDiagSCMatrix&).
        Implemented in terms of set_occupations(const RefDiagSCMatrix&,bool). */
    void set_occupations(const RefDiagSCMatrix& evals);
    virtual void set_occupations(const RefDiagSCMatrix& evals, bool can_change_multiplicity);    

    // scf things
    void init_vector();
    void done_vector();
    void reset_density();
    double new_density();
    double scf_energy();

    Ref<SCExtrapData> extrap_data();
    
    // gradient things
    void init_gradient();
    void done_gradient();

    RefSymmSCMatrix lagrangian();
    RefSymmSCMatrix gradient_density();

    // hessian things
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
