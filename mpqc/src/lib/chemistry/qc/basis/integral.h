//
// integral.h --- definition of the Integral class
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

#ifndef _chemistry_qc_basis_integral_h
#define _chemistry_qc_basis_integral_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/state/state.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>

class SymmetryOperation;
class RefPetiteList;
class RefSymmSCMatrix;
class ShellRotation;
class CartesianIter;
class RedundantCartesianIter;
class RedundantCartesianSubIter;
class SphericalTransformIter;
class SphericalTransform;
class PointBag_double;

SavableState_REF_fwddec(SCElementOp);
SavableState_REF_fwddec(SCElementOp3);

/** The Integral abstract class acts as a factory to provide objects that
compute one and two electron integrals.  */
class Integral : public SavableState {
#   define CLASSNAME Integral
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    /** Initialize the Integral object given a GaussianBasisSet for
        each center. */
    Integral(const RefGaussianBasisSet &b1,
             const RefGaussianBasisSet &b2,
             const RefGaussianBasisSet &b3,
             const RefGaussianBasisSet &b4);
    RefGaussianBasisSet bs1_;
    RefGaussianBasisSet bs2_;
    RefGaussianBasisSet bs3_;
    RefGaussianBasisSet bs4_;

    // the maximum number of bytes that should be used for
    // storing intermediates
    int storage_;
    int storage_used_;

    RefMessageGrp grp_;
  public:
    /// Restore the Integral object from the given StateIn object.
    Integral(StateIn&);
    /// Integral the Integral object from the given KeyVal object.
    Integral(const RefKeyVal&);
    
    void save_data_state(StateOut&);

    /// Sets the total amount of storage, in bytes, that is available.
    void set_storage(int i) { storage_=i; };
    /// Returns how much storage has been used.
    int storage_used() { return storage_used_; }
    /// Returns how much storage was not needed.
    int storage_unused();

    /** The specific integral classes use this to tell Integral
        how much memory they are using/freeing. */
    void adjust_storage(int s) { storage_used_ += s; }

    /// Return the PetiteList object.
    RefPetiteList petite_list();
    /// Return the PetiteList object for the given basis set.
    RefPetiteList petite_list(const RefGaussianBasisSet&);
    /** Return the ShellRotation object for a shell of the given angular
        momentum.  Pass nonzero to pure to do solid harmonics. */
    ShellRotation shell_rotation(int am, SymmetryOperation&, int pure=0);

    /// Set the basis set for each center.
    virtual void set_basis(const RefGaussianBasisSet &b1,
                           const RefGaussianBasisSet &b2 = 0,
                           const RefGaussianBasisSet &b3 = 0,
                           const RefGaussianBasisSet &b4 = 0);

    // /////////////////////////////////////////////////////////////////////
    // the following must be defined in the specific integral package

    /** Return a CartesianIter object.  The caller is responsible for
        freeing the object. */
    virtual CartesianIter * new_cartesian_iter(int) =0;
    /** Return a RedundantCartesianIter object.  The caller is responsible
        for freeing the object. */
    virtual RedundantCartesianIter * new_redundant_cartesian_iter(int) =0;
    /** Return a RedundantCartesianSubIter object.  The caller is
        responsible for freeing the object. */
    virtual RedundantCartesianSubIter *
                                 new_redundant_cartesian_sub_iter(int) =0;
    /** Return a SphericalTransformIter object.  The caller is
        responsible for freeing the object. */
    virtual SphericalTransformIter *
                  new_spherical_transform_iter(int l,
                                               int inv=0, int subl=-1) =0;
    /** Return a SphericalTransform object.  The pointer is only valid
        while this Integral object is valid. */
    virtual const SphericalTransform *
                  spherical_transform(int l,
                                      int inv=0, int subl=-1) =0;
    
    /// Return a OneBodyInt that computes the overlap.
    virtual RefOneBodyInt overlap() =0;
    
    /// Return a OneBodyInt that computes the kinetic energy.
    virtual RefOneBodyInt kinetic() =0;

    /** Return a OneBodyInt that computes the integrals for interactions
        with point charges. */
    virtual RefOneBodyInt point_charge(const RefPointChargeData&) =0;

    /** Return a OneBodyInt that computes the nuclear repulsion integrals.
        Charges from the atoms on the center one are used. */
    virtual RefOneBodyInt nuclear() = 0;

    /// Return a OneBodyInt that computes the core Hamiltonian integrals.
    virtual RefOneBodyInt hcore() = 0;

    /** Return a OneBodyInt that computes the electric field integrals
        dotted with a given vector. */
    virtual RefOneBodyInt efield_dot_vector(const RefEfieldDotVectorData&) =0;

    /// Return a OneBodyInt that computes dipole moment integrals.
    virtual RefOneBodyInt dipole(const RefDipoleData&) =0;

    /// Return a OneBodyDerivInt that computes overlap derivatives.
    virtual RefOneBodyDerivInt overlap_deriv() =0;
                                             
    /// Return a OneBodyDerivInt that computes kinetic energy derivatives.
    virtual RefOneBodyDerivInt kinetic_deriv() =0;
                                             
    /// Return a OneBodyDerivInt that computes nuclear repulsion derivatives.
    virtual RefOneBodyDerivInt nuclear_deriv() =0;
                                     
    /// Return a OneBodyDerivInt that computes core Hamiltonian derivatives.
    virtual RefOneBodyDerivInt hcore_deriv() =0;
                                             
    /// Return a TwoBodyInt that computes electron repulsion integrals.
    virtual RefTwoBodyInt electron_repulsion() =0;
    
    /// Return a TwoBodyDerivInt that computes electron repulsion derivatives.
    virtual RefTwoBodyDerivInt electron_repulsion_deriv() =0;

    /// Return the MessageGrp used by the integrals objects.
    RefMessageGrp messagegrp() { return grp_; }
};
SavableState_REF_dec(Integral);

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
