//
// wfn.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#ifndef _chemistry_qc_wfn_wfn_h
#define _chemistry_qc_wfn_wfn_h

#include <iostream>

#include <util/misc/compute.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <util/misc/xml.h>

using boost::property_tree::ptree;

namespace sc {

  /// @addtogroup ChemistryElectronicStructure
  /// @{

/** A Wavefunction is a MolecularEnergy that utilizies a GaussianBasisSet. */
class Wavefunction: public MolecularEnergy {

    RefSCDimension aodim_;
    RefSCDimension sodim_;
    Ref<SCMatrixKit> basiskit_;

    ResultRefSymmSCMatrix overlap_;
    ResultRefSymmSCMatrix hcore_;

  protected:
    ResultRefSCMatrix natural_orbitals_;
    ResultRefDiagSCMatrix natural_density_;
    Ref<GaussianBasisSet> gbs_;

  private:
    double * bs_values;
    double * bsg_values;

    Ref<Integral> integral_;

    Ref<GaussianBasisSet> atom_basis_;
    double * atom_basis_coef_;

    double lindep_tol_;
    OverlapOrthog::OrthogMethod orthog_method_;
    Ref<OverlapOrthog> orthog_;

    int print_nao_;
    int print_npa_;

    int dk_; // 0, 1, 2, or 3
    Ref<GaussianBasisSet> momentum_basis_;

    mutable double magnetic_moment_; //!< caches the value returned by magnetic_moment()

    void init_orthog();

    void set_up_charge_types(std::vector<int> &q_pc,
                             std::vector<int> &q_cd,
                             std::vector<int> &n_pc,
                             std::vector<int> &n_cd);

    double nuc_rep_pc_pc(const std::vector<int>&,const std::vector<int>&,bool);
    double nuc_rep_pc_cd(const std::vector<int>&,const std::vector<int>&);
    double nuc_rep_cd_cd(const std::vector<int>&,const std::vector<int>&,bool);
    void scale_atom_basis_coef();

    void nuc_rep_grad_pc_pc(double **grad,
                            const std::vector<int>&c1,
                            const std::vector<int>&c2,
                            bool uniq);
    void nuc_rep_grad_pc_cd(double **grad,
                            const std::vector<int>&c1,
                            const std::vector<int>&c2);
    void nuc_rep_grad_cd_cd(double **grad,
                            const std::vector<int>&c1,
                            const std::vector<int>&c2,
                            bool uniq);

    // Compute the dk relativistic core hamiltonian.
    RefSymmSCMatrix core_hamiltonian_dk(
      int dk,
      const Ref<GaussianBasisSet> &bas,
      const Ref<GaussianBasisSet> &pbas = 0);
    void core_hamiltonian_dk2_contrib(const RefSymmSCMatrix &h_pbas,
                                      const RefDiagSCMatrix &E,
                                      const RefDiagSCMatrix &K,
                                      const RefDiagSCMatrix &p2,
                                      const RefDiagSCMatrix &p2K2,
                                      const RefDiagSCMatrix &p2K2_inv,
                                      const RefSymmSCMatrix &AVA_pbas,
                                      const RefSymmSCMatrix &BpVpB_pbas);
    void core_hamiltonian_dk3_contrib(const RefSymmSCMatrix &h_pbas,
                                      const RefDiagSCMatrix &E,
                                      const RefDiagSCMatrix &B,
                                      const RefDiagSCMatrix &p2K2_inv,
                                      const RefSCMatrix &so_to_p,
                                      const RefSymmSCMatrix &pxVp);

    /// Wavefunction (reluctantly) supports calculations in finite electric fields in c1 symmetry
    /// general support is coming in the future.
    bool nonzero_efield_supported() const;

  protected:

    int debug_;

    double min_orthog_res();
    double max_orthog_res();

    void copy_orthog_info(const Ref<Wavefunction> &);

  public:
    Wavefunction(StateIn&);
    /** The KeyVal constructor. It accepts all keywords of MolecularEnergy and the following additional keywords:

        <dl>

        <dt><tt>basis</tt><dd> Specifies a GaussianBasisSet object.  There
        is no default.

        <dt><tt>integrals</tt><dd> Specifies an Integral object that
        computes the two electron integrals.  The default is a IntegralV3
        object.

        <dt><tt>orthog_method</tt><dd> This is a string that specifies the
        orthogonalization method to be used.  It can be one one canonical,
        gramschmidt, or symmetric.  The default is symmetric.

        <dt><tt>lindep_tol</tt><dd> The tolerance used to detect linearly
        dependent basis functions.  The precise meaning depends on the
        orthogonalization method.  The default value is 1e-8.

        <dt><tt>print_nao</tt><dd> This specifies a boolean value.  If true
        the natural atomic orbitals will be printed.  Not all wavefunction
        will be able to do this.  The default is false.

        <dt><tt>print_npa</tt><dd> This specifies a boolean value.  If true
        the natural population analysis will be printed.  Not all
        wavefunction will be able to do this.  The default is true if
        print_nao is true, otherwise it is false.

        <dt><tt>debug</tt><dd> This integer can be used to produce output
        for debugging.  The default is 0.

        </dl> */
    Wavefunction(const Ref<KeyVal>&);
    virtual ~Wavefunction();

    void save_data_state(StateOut&);
    virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

    double density(const SCVector3&);
    double density_gradient(const SCVector3&,double*);
    double natural_orbital(const SCVector3& r, int iorb);
    double natural_orbital_density(const SCVector3& r,
                                   int orb, double* orbval = 0);
    double orbital(const SCVector3& r, int iorb, const RefSCMatrix& orbs);
    // same as above, except computes all orbital AND orbs must be transposed (i.e. mo rows, ao columns)
    void orbitals(const SCVector3& r, const RefSCMatrix& orbs, RefSCVector& values);
    // same as above except it evaluates on a set of values, and is static
    static void orbitals(const Ref<OrbitalSpace>& orbs,
                         const std::vector<SCVector3>& r,
                         std::vector<double>& values);

    double orbital_density(const SCVector3& r,
                           int iorb,
                           const RefSCMatrix& orbs,
                           double* orbval = 0);

    /// Returns the total charge of the system.
    double total_charge() const;
    /// Returns the number of electrons.
    virtual int nelectron() = 0;
    /// Computes the S (or J) magnetic moment
    /// of the target state(s), in units of \f$ \hbar/2 \f$.
    /// Can be evaluated from density and overlap, as;
    /// \code
    ///   (this->alpha_density() * this-> overlap()).trace() -
    ///   (this->beta_density() * this-> overlap()).trace()
    /// \endcode
    /// but derived Wavefunction may have this value as user input.
    /// @return the magnetic moment
    virtual double magnetic_moment() const;

    /// Returns the SO density.
    virtual RefSymmSCMatrix density() = 0;
    /// Returns the AO density.
    virtual RefSymmSCMatrix ao_density();
    /// Returns the natural orbitals, in SO basis
    virtual RefSCMatrix natural_orbitals();
    /// Returns the natural density (a diagonal matrix).
    virtual RefDiagSCMatrix natural_density();

    /// Return 1 if the magnetic moment != 0
    int spin_polarized() { return magnetic_moment() != 0.0; }

    /// Returns the level the of the Douglas-Kroll approximation.
    int dk() const { return dk_; }

    /// Return alpha electron densities in the SO basis.
    virtual RefSymmSCMatrix alpha_density();
    /// Return beta electron densities in the SO basis.
    virtual RefSymmSCMatrix beta_density();
    /// Return alpha electron densities in the AO basis.
    virtual RefSymmSCMatrix alpha_ao_density();
    /// Return beta electron densities in the AO basis.
    virtual RefSymmSCMatrix beta_ao_density();

    /// returns the ao to nao transformation matrix
    virtual RefSCMatrix nao(double *atom_charges=0);

    /// Returns the SO overlap matrix.
    virtual RefSymmSCMatrix overlap();
    /** Returns the SO core Hamiltonian in the given basis and momentum
        basis.  The momentum basis is not needed if no Douglas-Kroll
        correction is being performed. */
    virtual RefSymmSCMatrix core_hamiltonian_for_basis(
      const Ref<GaussianBasisSet> &bas,
      const Ref<GaussianBasisSet> &pbas = 0);
    /// Returns the SO core Hamiltonian.
    virtual RefSymmSCMatrix core_hamiltonian();
    // Compute the non-relativistic core hamiltonian.
    RefSymmSCMatrix core_hamiltonian_nr(
      const Ref<GaussianBasisSet> &bas);


    /** Returns the nuclear repulsion energy.  This must be used instead of
        Molecule::nuclear_repulsion_energy() since there may be contributions from
        diffuse atomic charges or an external electric field (\sa MolecularEnergy::electric_field). */
    virtual double nuclear_repulsion_energy();
    /** Computes the nuclear repulsion gradient.  This must be used instead
        of Molecule::nuclear_repulsion_1der() since there may be contributions from diffuse
        atomic charges.  The gradient, g, is zeroed and set to x_0, y_0,
        z_0, x_1, ... . */
    void nuclear_repulsion_energy_gradient(double *g);
    /** Computes the nuclear repulsion gradient.  This must be used instead
        of Molecule::nuclear_repulsion_1der() since there may be contributions from diffuse
        atomic charges.  The gradient, g, is first zeroed.  Its dimensions
        are g[natom][3]. */
    virtual void nuclear_repulsion_energy_gradient(double **g);

    /// Atomic orbital dimension.
    RefSCDimension ao_dimension();
    /// Symmetry adapted orbital dimension.
    RefSCDimension so_dimension();
    /// Orthogonalized symmetry adapted orbital dimension.
    RefSCDimension oso_dimension();
    /// Matrix kit for AO, SO, orthogonalized SO, and MO dimensioned matrices.
    Ref<SCMatrixKit> basis_matrixkit();
    /// Returns the Molecule.
    Ref<Molecule> molecule() const;
    /// Returns the basis set.
    Ref<GaussianBasisSet> basis() const;
    /// Returns the basis used for p^2 in the DK correction.
    Ref<GaussianBasisSet> momentum_basis() const;
    /// Returns the basis set describing the nuclear charge distributions
    Ref<GaussianBasisSet> atom_basis() const;
    /** Returns the coefficients of the nuclear charge distribution basis
     * functions. */
    const double *atom_basis_coef() const;
    /// Returns the integral evaluator.
    Ref<Integral> integral();

    // override symmetry_changed from MolecularEnergy
    void symmetry_changed();

    /** Returns a matrix which does the default transform from SO's to
        orthogonal SO's.  This could be either the symmetric or canonical
        orthogonalization matrix.  The row dimension is SO and the column
        dimension is ortho SO.  An operator \f$O\f$ in the ortho SO basis is
        given by \f$X O X^T\f$ where \f$X\f$ is the return value of this
        function. */
    RefSCMatrix so_to_orthog_so();

    /** Returns the inverse of the transformation returned by so_to_orthog_so.
     */
    RefSCMatrix so_to_orthog_so_inverse();

    /// Returns the orthogonalization method
    OverlapOrthog::OrthogMethod orthog_method() const;
    /** (Re)Sets the orthogonalization method and makes this obsolete.
        Must be virtual to allow proper obsoletion of the object.
     */
    virtual void set_orthog_method(const OverlapOrthog::OrthogMethod&);

    /// Returns the tolerance for linear dependencies.
    double lindep_tol() const;
    /// Re(Sets) the tolerance for linear dependencies.
    void set_lindep_tol(double);

    void obsolete();

    void print(std::ostream& = ExEnv::out0()) const;


    /** output orbitals to some files to facilitate plotting, with the help of the WriteOrbital class.
     */
    void writeorbitals();
};

/// @}
// end of addtogroup ChemistryElectronicStructure

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
