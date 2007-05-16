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

#include <stddef.h>

#include <util/state/state.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/intparams.h>

namespace sc {

class SymmetryOperation;
class RefSymmSCMatrix;
class ShellRotation;
class CartesianIter;
class RedundantCartesianIter;
class RedundantCartesianSubIter;
class SphericalTransformIter;
class SphericalTransform;
class PointBag_double;
class PetiteList;

/** The Integral abstract class acts as a factory to provide objects that
compute one and two electron integrals.  */
class Integral : public SavableState {
  protected:
    /** Initialize the Integral object given a GaussianBasisSet for
        each center. */
    Integral(const Ref<GaussianBasisSet> &b1,
             const Ref<GaussianBasisSet> &b2,
             const Ref<GaussianBasisSet> &b3,
             const Ref<GaussianBasisSet> &b4);
    Ref<GaussianBasisSet> bs1_;
    Ref<GaussianBasisSet> bs2_;
    Ref<GaussianBasisSet> bs3_;
    Ref<GaussianBasisSet> bs4_;

    // the maximum number of bytes that should be used for
    // storing intermediates
    size_t storage_;
    size_t storage_used_;

    Ref<MessageGrp> grp_;
  public:
    /// Restore the Integral object from the given StateIn object.
    Integral(StateIn&);
    /// Construct the Integral object from the given KeyVal object.
    Integral(const Ref<KeyVal>&);

    virtual ~Integral();
    
    void save_data_state(StateOut&);

    /** Create an integral factory.  This routine looks for a -integral
        argument, then the environmental variable INTEGRAL.
        The argument to -integral should
        be either string for a ParsedKeyVal constructor or a classname.
        This factory is not guaranteed to have its storage and basis
        sets set up properly, hence set_basis and set_storage
        need to be called on it. */
    static Integral* initial_integral(int &argc, char **argv);
    /// Specifies a new default Integral factory
    static void set_default_integral(const Ref<Integral>&);
    /// Returns the default Integral factory
    static Integral* get_default_integral();
    /// Clones the given Integral factory. The new factory may need to have set_basis and set_storage to be called on it.
    virtual Integral* clone() =0;

    /** Returns nonzero if this and the given Integral object have the same
        integral ordering, normalization conventions, etc.  */
    virtual int equiv(const Ref<Integral> &);

    /// Sets the total amount of storage, in bytes, that is available.
    virtual void set_storage(size_t i) { storage_=i; };
    /// Returns how much storage has been used.
    size_t storage_used() { return storage_used_; }
    /// Returns how much storage was not needed.
    size_t storage_unused();
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for electron repulsion integrals. */
    virtual size_t storage_required_eri(const Ref<GaussianBasisSet> &b1,
					const Ref<GaussianBasisSet> &b2 = 0,
					const Ref<GaussianBasisSet> &b3 = 0,
					const Ref<GaussianBasisSet> &b4 = 0);
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for linear R12 integrals. */
    virtual size_t storage_required_grt(const Ref<GaussianBasisSet> &b1,
					const Ref<GaussianBasisSet> &b2 = 0,
					const Ref<GaussianBasisSet> &b3 = 0,
					const Ref<GaussianBasisSet> &b4 = 0);
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for G12 integrals. */
    virtual size_t storage_required_g12(const Ref<GaussianBasisSet> &b1,
					const Ref<GaussianBasisSet> &b2 = 0,
					const Ref<GaussianBasisSet> &b3 = 0,
					const Ref<GaussianBasisSet> &b4 = 0);
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for G12 integrals. */
    virtual size_t storage_required_g12nc(const Ref<GaussianBasisSet> &b1,
					  const Ref<GaussianBasisSet> &b2 = 0,
					  const Ref<GaussianBasisSet> &b3 = 0,
					  const Ref<GaussianBasisSet> &b4 = 0);
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for general G12 integrals. */
    virtual size_t storage_required_geng12(const Ref<GaussianBasisSet> &b1,
					   const Ref<GaussianBasisSet> &b2 = 0,
					   const Ref<GaussianBasisSet> &b3 = 0,
					   const Ref<GaussianBasisSet> &b4 = 0);
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for derivative electron repulsion integrals. */
    virtual size_t storage_required_eri_deriv(const Ref<GaussianBasisSet> &b1,
					      const Ref<GaussianBasisSet> &b2 = 0,
					      const Ref<GaussianBasisSet> &b3 = 0,
					      const Ref<GaussianBasisSet> &b4 = 0);

    /** The specific integral classes use this to tell Integral
        how much memory they are using/freeing. */
    void adjust_storage(ptrdiff_t s) { storage_used_ += s; }

    /// Return the PetiteList object.
    Ref<PetiteList> petite_list();
    /// Return the PetiteList object for the given basis set.
    Ref<PetiteList> petite_list(const Ref<GaussianBasisSet>&);
    /** Return the ShellRotation object for a shell of the given angular
        momentum.  Pass nonzero to pure to do solid harmonics. */
    ShellRotation shell_rotation(int am, SymmetryOperation&, int pure=0);

    /// Set the basis set for each center.
    virtual void set_basis(const Ref<GaussianBasisSet> &b1,
                           const Ref<GaussianBasisSet> &b2 = 0,
                           const Ref<GaussianBasisSet> &b3 = 0,
                           const Ref<GaussianBasisSet> &b4 = 0);

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
    virtual RedundantCartesianSubIter*
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
    virtual Ref<OneBodyInt> overlap() =0;
    
    /// Return a OneBodyInt that computes the kinetic energy.
    virtual Ref<OneBodyInt> kinetic() =0;

    /** Return a OneBodyInt that computes the integrals for interactions
        with point charges. */
    virtual Ref<OneBodyInt> point_charge(const Ref<PointChargeData>&) =0;

    /** Return a OneBodyInt that computes the integrals for interactions
        with point charges. */
    virtual Ref<OneBodyOneCenterInt> point_charge1(const Ref<PointChargeData>&);

    /** Return a OneBodyInt that computes the nuclear repulsion integrals.
        Charges from the atoms on center one are used.  If center two is
        not identical to center one, then the charges on center two are
        included as well.  */
    virtual Ref<OneBodyInt> nuclear() = 0;

    /** Return a OneBodyInt that computes $\bar{p}\cdotp V\bar{p}$, where
        $V$ is the nuclear potential. */
    virtual Ref<OneBodyInt> p_dot_nuclear_p();

    /** Return a OneBodyInt that computes $\bar{p}\times V\bar{p}$, where
        $V$ is the nuclear potential. This is different than most other
        one body integrals, in that each entry in the integral buffer
        is a vector of three integrals. */
    virtual Ref<OneBodyInt> p_cross_nuclear_p();

    /// Return a OneBodyInt that computes the core Hamiltonian integrals.
    virtual Ref<OneBodyInt> hcore() = 0;

    /** Return a OneBodyInt that computes the electric field integrals
        dotted with a given vector. */
    virtual Ref<OneBodyInt> efield_dot_vector(const Ref<EfieldDotVectorData>&) =0;

    /** Return a OneBodyInt that computes electric dipole moment integrals.
	The canonical order of integrals in a set is x, y, z. */
    virtual Ref<OneBodyInt> dipole(const Ref<DipoleData>&) =0;

    /** Return a OneBodyInt that computes electric quadrupole moment integrals.
	The canonical order of integrals in a set is x^2, xy, xz, y^2, yz, z^2. */
    virtual Ref<OneBodyInt> quadrupole(const Ref<DipoleData>&) =0;

    /// Return a OneBodyDerivInt that computes overlap derivatives.
    virtual Ref<OneBodyDerivInt> overlap_deriv() =0;
                                             
    /// Return a OneBodyDerivInt that computes kinetic energy derivatives.
    virtual Ref<OneBodyDerivInt> kinetic_deriv() =0;
                                             
    /// Return a OneBodyDerivInt that computes nuclear repulsion derivatives.
    virtual Ref<OneBodyDerivInt> nuclear_deriv() =0;
                                     
    /// Return a OneBodyDerivInt that computes core Hamiltonian derivatives.
    virtual Ref<OneBodyDerivInt> hcore_deriv() =0;

    /** Return a TwoBodyThreeCenterInt that computes electron repulsion
        integrals.  If this is not re-implemented it will throw. */
    virtual Ref<TwoBodyThreeCenterInt> electron_repulsion3();

    /** Return a TwoBodyThreeCenterInt that computes electron repulsion
        integrals.  If this is not re-implemented it will throw. */
    virtual Ref<TwoBodyThreeCenterDerivInt> electron_repulsion3_deriv();

    /** Return a TwoBodyTwoCenterInt that computes electron repulsion
        integrals. If this is not re-implemented it will throw. */
    virtual Ref<TwoBodyTwoCenterInt> electron_repulsion2();

    /** Return a TwoBodyTwoCenterInt that computes electron repulsion
        integrals. If this is not re-implemented it will throw. */
    virtual Ref<TwoBodyTwoCenterDerivInt> electron_repulsion2_deriv();

    /// Return a TwoBodyInt that computes electron repulsion integrals.
    virtual Ref<TwoBodyInt> electron_repulsion();

    /// Return a TwoBodyDerivInt that computes electron repulsion derivatives.
    virtual Ref<TwoBodyDerivInt> electron_repulsion_deriv();

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to linear R12 methods.  According to the convention in the
        literature, "g" stands for electron repulsion integral, "r" for the
        integral of r12 operator, and "t" for the commutator
        integrals. Implementation for this kind of TwoBodyInt is
        optional. */
    virtual Ref<TwoBodyInt> grt();
    
    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use Gaussian geminals.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> g12(const Ref<IntParamsG12>&);

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use Gaussian geminals.
	This particular implementation does not produce commutator integrals
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> g12nc(const Ref<IntParamsG12>&);

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use general Gaussian geminals (i.e. exp(-a*(r1+r2)-g*r12)).
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> geng12(const Ref<IntParamsGenG12>&);
    
    /// Return the MessageGrp used by the integrals objects.
    Ref<MessageGrp> messagegrp() { return grp_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
