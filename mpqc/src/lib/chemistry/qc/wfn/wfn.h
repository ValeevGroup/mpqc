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

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <util/misc/compute.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>

/** A Wavefunction is a MolecularEnergy that utilizies a GaussianBasisSet. */
class Wavefunction: public MolecularEnergy {
#   define CLASSNAME Wavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefSCDimension aodim_;
    RefSCDimension sodim_;
    RefSCDimension osodim_;
    RefSCMatrixKit basiskit_;

    int have_overlap_eig_;
    RefSCMatrix overlap_eigvec_;
    RefDiagSCMatrix overlap_isqrt_eigval_;
    RefDiagSCMatrix overlap_sqrt_eigval_;
    
    ResultRefSymmSCMatrix overlap_;
    ResultRefSymmSCMatrix hcore_;
    ResultRefSCMatrix natural_orbitals_;
    ResultRefDiagSCMatrix natural_density_;

    double * bs_values;
    double * bsg_values;

    RefGaussianBasisSet gbs_;
    RefIntegral integral_;

    // The tolerance for lambda(max)/lambda(min) for linearly
    // independent basis functions
    double lindep_tol_;

    // Whether or not to symmetrically orthogonalize
    int symm_orthog_;

    int print_nao_;
    int print_npa_;

    void compute_overlap_eig();

  protected:

    int debug_;
    
  public:
    Wavefunction(StateIn&);
    /** @memo The KeyVal constructor.

        \begin{description}

        \item[basis] Specifies a GaussianBasisSet object.  There is no
        default.

        \item[integral] Specifies an Integral object that computes the two
        electron integrals.  The default is a IntegralV3 object.

        \item[symm_orthog] If true, symmetric orthogonalization is used;
        otherwise canonical orthogonalization is used.  The default is true.

        \item[print_nao] This specifies a boolean value.  If true the
        natural atomic orbitals will be printed.  Not all wavefunction will
        be able to do this.  The default is false.

        \item[print_npa] This specifies a boolean value.  If true the
        natural population analysis will be printed.  Not all wavefunction
        will be able to do this.  The default is true if print_nao is true,
        otherwise it is false.

        \item[debug] This integer can be used to produce output for
        debugging.  The default is 0.

        \end{description}
     */
    Wavefunction(const RefKeyVal&);
    virtual ~Wavefunction();

    void save_data_state(StateOut&);

    double density(const SCVector3&);
    double density_gradient(const SCVector3&,double*);
    double natural_orbital(const SCVector3& r, int iorb);
    double natural_orbital_density(const SCVector3& r,
                                   int orb, double* orbval = 0);
    double orbital(const SCVector3& r, int iorb, const RefSCMatrix& orbs);

    double orbital_density(const SCVector3& r,
                           int iorb,
                           const RefSCMatrix& orbs,
                           double* orbval = 0);

    /// Returns the charge.
    double charge();
    /// Returns the number of electrons.
    virtual int nelectron() = 0;

    /// Returns the SO density.
    virtual RefSymmSCMatrix density() = 0;
    /// Returns the AO density.
    virtual RefSymmSCMatrix ao_density();
    /// Returns the natural orbitals.
    virtual RefSCMatrix natural_orbitals();
    /// Returns the natural density (a diagonal matrix).
    virtual RefDiagSCMatrix natural_density();

    /// Return 1 if the alpha density is not equal to the beta density.
    virtual int spin_polarized() = 0;

    /// Return alpha electron densities in the SO basis.
    virtual RefSymmSCMatrix alpha_density();
    /// Return beta electron densities in the SO basis.
    virtual RefSymmSCMatrix beta_density();
    /// Return alpha electron densities in the AO basis.
    virtual RefSymmSCMatrix alpha_ao_density();
    /// Return beta electron densities in the AO basis.
    virtual RefSymmSCMatrix beta_ao_density();

    /// returns the ao to nao transformation matrix
    virtual RefSCMatrix nao();

    /// Returns the SO overlap matrix.
    virtual RefSymmSCMatrix overlap();
    /// Returns the SO core Hamiltonian.
    virtual RefSymmSCMatrix core_hamiltonian();

    /// Atomic orbital dimension.
    RefSCDimension ao_dimension();
    /// Symmetry adapted orbital dimension.
    RefSCDimension so_dimension();
    /// Orthogonalized symmetry adapted orbital dimension.
    RefSCDimension oso_dimension();
    /// Matrix kit for AO, SO, orthogonalized SO, and MO dimensioned matrices.
    RefSCMatrixKit basis_matrixkit();
    /// Returns the basis set.
    RefGaussianBasisSet basis() const;
    /// Returns the integral evaluator.
    RefIntegral integral();

    // override symmetry_changed from MolecularEnergy
    void symmetry_changed();

    /** Returns a matrix which does the default transform from SO's to
        orthogonal SO's.  This could be either the symmetric or canonical
        orthogonalization matrix.  The row dimension is SO and the column
        dimension is ortho SO.  An operator $O$ in the ortho SO basis is
        given by $X O X^T$ where $X$ is the return value of this
        function. */
    RefSCMatrix so_to_orthog_so();

    /** Returns the inverse of the transformation returned by so_to_orthog_so.
     */
    RefSCMatrix so_to_orthog_so_inverse();

    void print(ostream& = cout) const;
};
SavableState_REF_dec(Wavefunction);

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
