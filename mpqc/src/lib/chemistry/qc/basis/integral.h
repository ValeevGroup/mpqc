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
class PointBag_double;

SavableState_REF_fwddec(SCElementOp);
SavableState_REF_fwddec(SCElementOp3);

// some useful things to have that depend on the underlying integrals
// package

class Integral : public SavableState {
#   define CLASSNAME Integral
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
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

    RefMessageGrp grp_;

  public:
    Integral(StateIn&);
    Integral(const RefKeyVal&);
    
    void save_data_state(StateOut&);

    void set_storage(int i) { storage_=i; };

    RefPetiteList petite_list();
    RefPetiteList petite_list(const RefGaussianBasisSet&);
    ShellRotation shell_rotation(int am, SymmetryOperation&, int pure=0);

    virtual void set_basis(const RefGaussianBasisSet &b1,
                           const RefGaussianBasisSet &b2 = 0,
                           const RefGaussianBasisSet &b3 = 0,
                           const RefGaussianBasisSet &b4 = 0);

    ///////////////////////////////////////////////////////////////////////
    // the following must be defined in the specific integral package

    virtual CartesianIter * new_cartesian_iter(int) =0;
    virtual RedundantCartesianIter * new_redundant_cartesian_iter(int) =0;
    virtual RedundantCartesianSubIter *
                                 new_redundant_cartesian_sub_iter(int) =0;
    virtual SphericalTransformIter *
                              new_spherical_transform_iter(int, int=0) =0;
    
    virtual RefOneBodyInt overlap() =0;
    
    virtual RefOneBodyInt kinetic() =0;

    virtual RefOneBodyInt point_charge(const RefPointChargeData&) =0;

    // charges from the atom on the first center are used
    virtual RefOneBodyInt nuclear() = 0;

    virtual RefOneBodyInt hcore() = 0;

    virtual RefOneBodyInt efield_dot_vector(const RefEfieldDotVectorData&) =0;

    virtual RefOneBodyInt dipole(const RefDipoleData&) =0;

    virtual RefOneBodyDerivInt overlap_deriv() =0;
                                             
    virtual RefOneBodyDerivInt kinetic_deriv() =0;
                                             
    virtual RefOneBodyDerivInt nuclear_deriv() =0;
                                     
    virtual RefOneBodyDerivInt hcore_deriv() =0;
                                             
    virtual RefTwoBodyInt electron_repulsion() =0;
    
    virtual RefTwoBodyDerivInt electron_repulsion_deriv() =0;
};
SavableState_REF_dec(Integral);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
