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

#include <stddef.h>

#include <util/state/state.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/intparams.h>
#include <chemistry/qc/basis/operator.h>

#define INTEGRALTHROWDEFAULT { throw ProgrammingError("not implemented in this particular factory", __FILE__, __LINE__); }

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

/// @addtogroup ChemistryBasisIntegral
/// @{

/// returns the type of the evaluator for evaluating this set of two-body integrals
template <int NumCenters> struct TwoBodyIntEvalType;
template <> struct TwoBodyIntEvalType<4> {
  typedef TwoBodyInt value;
};
template <> struct TwoBodyIntEvalType<3> {
  typedef TwoBodyThreeCenterInt value;
};
template <> struct TwoBodyIntEvalType<2> {
  typedef TwoBodyTwoCenterInt value;
};
// forward declarations of auxiliary type functions
namespace detail {
  template <int NumCenters> struct ERIEvalCreator;
  template <int NumCenters> struct R12EvalCreator;
  template <int NumCenters> struct G12EvalCreator;
  template <int NumCenters> struct G12NCEvalCreator;
  template <int NumCenters> struct G12DKHEvalCreator;
  template <int NumCenters> struct R120G12EvalCreator;
  template <int NumCenters> struct R12m1G12EvalCreator;
  template <int NumCenters> struct G12T1G12EvalCreator;
  template <int NumCenters> struct DeltaFunctionEvalCreator;
  template <int NumCenters, TwoBodyOperSet::type Type> struct EvalCreator;
};


/** The Integral abstract class acts as a factory to provide objects that compute one and two electron integrals.

  */
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

    typedef enum {
	MPQCSolidHarmonicsOrdering,
	CCASolidHarmonicsOrdering
    } SolidHarmonicsOrdering;
    SolidHarmonicsOrdering sharmorder_;

    // the maximum number of bytes that should be used for
    // storing intermediates
    size_t storage_;
    size_t storage_used_;
    // in absence of atomic<> weapons use this to protect storage_ and storage_used_
    Ref<ThreadLock> tlock_;

    Ref<MessageGrp> grp_;

  private:
    /**
       CCA standard specifies the new ordering (-l ... +l).
       MPQCSolidHarmonicsOrdering was the old ordering (0, +1, -1, +2, -2, etc.) and can no longer be used.
    */
    static const SolidHarmonicsOrdering default_sharmorder_ = CCASolidHarmonicsOrdering;

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
    /// Describes the ordering of the cartesian functions in a shell
    typedef enum {
    IntV3CartesianOrdering,
    CCACartesianOrdering,
    GAMESSCartesianOrdering,
    } CartesianOrdering;
    /// returns the ordering used by this factory
    virtual CartesianOrdering cartesian_ordering() const =0;

    /// Sets the total amount of storage, in bytes, that is available.
    virtual void set_storage(size_t i);
    /// Returns how much storage has been used.
    size_t storage_used() const { return storage_used_; }
    /// Returns how much storage was not needed.
    size_t storage_unused() const;

    /**
     * Reports the approximate amount of memory required, in bytes, to create an evaluator
     * for \c opertype
     *
     * @param opertype the operator type, TwoBodyOper::type
     * @param tbinttype the integral type, TwoBodyIntShape::value
     * @param deriv_level derivative level
     * @param b1 basis set on center 1
     * @param b2 basis set on center 2
     * @param b3 basis set on center 3
     * @param b4 basis set on center 4
     * @return number of bytes needed to create an evaluator of the specified type
     */
    virtual size_t storage_required(TwoBodyOper::type opertype,
                                    TwoBodyIntShape::value tbinttype,
                                    size_t deriv_level = 0,
                                    const Ref<GaussianBasisSet> &b1 = 0,
                                    const Ref<GaussianBasisSet> &b2 = 0,
                                    const Ref<GaussianBasisSet> &b3 = 0,
                                    const Ref<GaussianBasisSet> &b4 = 0);

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
      evaluator for G12NC integrals. */
    virtual size_t storage_required_g12nc(const Ref<GaussianBasisSet> &b1,
					  const Ref<GaussianBasisSet> &b2 = 0,
					  const Ref<GaussianBasisSet> &b3 = 0,
					  const Ref<GaussianBasisSet> &b4 = 0);
  /** Returns how much storage will be needed to initialize a two-body integrals
      evaluator for G12DKH integrals. */
    virtual size_t storage_required_g12dkh(const Ref<GaussianBasisSet> &b1,
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
    void adjust_storage(ptrdiff_t s);

    /// Return the PetiteList object.
    Ref<PetiteList> petite_list();
    /// Return the PetiteList object for the given basis set.
    Ref<PetiteList> petite_list(const Ref<GaussianBasisSet>&);
    /** Return the ShellRotation object for a shell of the given angular
        momentum.  Pass nonzero to pure to do solid harmonics. */
    ShellRotation shell_rotation(int am, SymmetryOperation&, int pure=0);

    /// retrieves basis for center 1
    const Ref<GaussianBasisSet>& basis1() const { return bs1_; }
    /// retrieves basis for center 2
    const Ref<GaussianBasisSet>& basis2() const { return bs2_; }
    /// retrieves basis for center 3
    const Ref<GaussianBasisSet>& basis3() const { return bs3_; }
    /// retrieves basis for center 4
    const Ref<GaussianBasisSet>& basis4() const { return bs4_; }
    /** Set the basis set for each center.
     * @param[in] b1 basis set on center 1; there is no default
     * @param[in] b2 basis set on center 2; if null, will use b1
     * @param[in] b3 basis set on center 3; if null, will use b2
     * @param[in] b4 basis set on center 4; if null, will use b3
     */
    virtual void set_basis(const Ref<GaussianBasisSet> &b1,
                           const Ref<GaussianBasisSet> &b2 = 0,
                           const Ref<GaussianBasisSet> &b3 = 0,
                           const Ref<GaussianBasisSet> &b4 = 0);

    /// Return the MessageGrp used by the integrals objects.
    Ref<MessageGrp> messagegrp() { return grp_; }

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
    /** Return a SphericalTransformIter object.  This factory must have been
        initialized with a basis set whose maximum angular momentum is greater than
        or equal to l. The caller is
        responsible for freeing the object. */
    virtual SphericalTransformIter *
                  new_spherical_transform_iter(int l,
                                               int inv=0, int subl=-1) =0;
    /** Return a SphericalTransform object.  This factory must have been
        initialized with a basis set whose maximum angular momentum is greater than
        or equal to l. The pointer is only valid
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

    /** Return a OneBodyInt that computes \f$\bar{p}\cdot V\bar{p}\f$, where
        \f$V\f$ is the nuclear potential. */
    virtual Ref<OneBodyInt> p_dot_nuclear_p();

    /** Return a OneBodyInt that computes \f$\bar{p}\times V\bar{p}\f$, where
        \f$V\f$ is the nuclear potential. This is different than most other
        one body integrals, in that each entry in the integral buffer
        is a vector of three integrals. */
    virtual Ref<OneBodyInt> p_cross_nuclear_p();

    /// Return a OneBodyInt that computes \f$p^4 = (\bar{p} \cdot \bar{p})^2\f$.
    virtual Ref<OneBodyInt> p4() =0;

    /// Return a OneBodyInt that computes the core Hamiltonian integrals.
    virtual Ref<OneBodyInt> hcore() = 0;

    /** Return a OneBodyInt that computes the electric field integrals at specified point.
        The canonical order of integrals in a set is x, y, z (i.e. Ex, Ey, Ey).

        @param O IntParamsOrigin object that specifies the point where the electric field is computed;
                 there is no default.

        \sa efield_dot_vector()
      */
    virtual Ref<OneBodyInt> efield(const Ref<IntParamsOrigin>& O) =0;

    /** Return a OneBodyInt that computes the electric field integrals at a given position
        dotted with a given vector. \sa efield()
      */
    virtual Ref<OneBodyInt> efield_dot_vector(const Ref<EfieldDotVectorData>&) =0;

    /** Return a OneBodyInt that computes the electric field gradient integrals at specified point.
        The canonical order of integrals in the 6-element sequence
        is d Ex / dx, d Ex / dy, d Ex / dz, d Ey / dy, d Ey / dz, d Ez / dz,

        @param O IntParamsOrigin object that specifies the point where the electric field gradient is computed;
                 there is no default.

        \note only 6 elements are unique since d Ei / d j = d Ej / d i

      */
    virtual Ref<OneBodyInt> efield_gradient(const Ref<IntParamsOrigin>& O);

    /** Return a OneBodyInt that computes electric dipole moment integrals, i.e. integrals
        of the \f$e (\mathbf{r}-\mathbf{O})\f$ operator. The canonical order of integrals in a set is x, y, z.

        @param O IntParamsOrigin object that specifies the origin of the multipole expansion;
                 the default is to use the origin of the coordinate system.

        \note Multiply by -1 to obtain electronic electric quadrupole integrals.
      */
    virtual Ref<OneBodyInt> dipole(const Ref<IntParamsOrigin>& O = 0) =0;

    /** Return a OneBodyInt that computes electric quadrupole moment integrals,
        i.e. integrals of the \f$e (\mathbf{r}-\mathbf{O}) \otimes (\mathbf{r}-\mathbf{O})\f$ operator.
	    The canonical order of integrals in a set is x^2, xy, xz, y^2, yz, z^2.

	    @param O IntParamsOrigin object that specifies the origin of the multipole expansion;
                 the default is to use the origin of the coordinate system

	    \note These are not traceless quadrupole integrals!!
	    \note Multiply by -1 to obtain electronic electric quadrupole integrals.

	  */
    virtual Ref<OneBodyInt> quadrupole(const Ref<IntParamsOrigin>& O  = 0) =0;

    /// Return a OneBodyDerivInt that computes overlap derivatives.
    virtual Ref<OneBodyDerivInt> overlap_deriv() =0;

    /// Return a OneBodyDerivInt that computes kinetic energy derivatives.
    virtual Ref<OneBodyDerivInt> kinetic_deriv() =0;

    /// Return a OneBodyDerivInt that computes nuclear repulsion derivatives.
    virtual Ref<OneBodyDerivInt> nuclear_deriv() =0;

    /// Return a OneBodyDerivInt that computes core Hamiltonian derivatives.
    virtual Ref<OneBodyDerivInt> hcore_deriv() =0;

    /** Return a TwoBodyThreeCenterInt that computes electron repulsion
        integrals. Electron 1 corresponds to centers 1 and 2, electron 2
        corresponds to center 3. If this is not re-implemented it will throw.

        @deprecated Use sc::Integral::coulomb<3>() instead.
        */
    DEPRECATED virtual Ref<TwoBodyThreeCenterInt> electron_repulsion3();

    /** Return a TwoBodyThreeCenterInt that computes electron repulsion
        integrals.  If this is not re-implemented it will throw. \sa electron_repulsion3() */
    virtual Ref<TwoBodyThreeCenterDerivInt> electron_repulsion3_deriv();

    /** Return a TwoBodyTwoCenterInt that computes electron repulsion
        integrals. If this is not re-implemented it will throw.

        @deprecated Use sc::Integral::coulomb<2>() instead.
        */
    DEPRECATED virtual Ref<TwoBodyTwoCenterInt> electron_repulsion2();

    /** Return a TwoBodyTwoCenterInt that computes electron repulsion
        integrals. If this is not re-implemented it will throw. */
    virtual Ref<TwoBodyTwoCenterDerivInt> electron_repulsion2_deriv();

    /** Return a TwoBodyInt that computes electron repulsion integrals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrERI.

        @deprecated Use sc::Integral::coulomb<4>() instead.
        */
    DEPRECATED virtual Ref<TwoBodyInt> electron_repulsion();

    /// Return a TwoBodyDerivInt that computes electron repulsion derivatives.
    virtual Ref<TwoBodyDerivInt> electron_repulsion_deriv();

    /**
     * Creates an evaluator for \c opertype
     *
     * @param opertype the operator type, TwoBodyOper::type
     * @param tbinttype the integral type, TwoBodyIntShape::value
     * @param deriv_level derivative level
     * @param b1 basis set on center 1
     * @param b2 basis set on center 2
     * @param b3 basis set on center 3
     * @param b4 basis set on center 4
     * @return number of bytes needed to create an evaluator of the specified type
     */
    virtual Ref<TwoBodyIntEval> make_eval(TwoBodyOper::type opertype,
                                          TwoBodyIntShape::value tbinttype,
                                          size_t deriv_level = 0,
                                          const Ref<GaussianBasisSet> &b1 = 0,
                                          const Ref<GaussianBasisSet> &b2 = 0,
                                          const Ref<GaussianBasisSet> &b3 = 0,
                                          const Ref<GaussianBasisSet> &b4 = 0);

    /** Return the evaluator of two-body integrals with Coulomb kernel:
        \f$
          r_{12}^{-1},
        \f$
        The evaluator will produce a set of integrals described by TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::ERI>.

        @tparam NumCenters specifies the number of centers that carry basis functions.
        Valid values are 4, 3, and 2.

        @note Implementation of this function is optional. The default implementation will throw FeatureNotImplemented . It is implemented in sc::IntegralV3 and sc::IntegralLibint2 .
      */
    template <int NumCenters>
    Ref< typename TwoBodyIntEvalType<NumCenters>::value > coulomb() {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::ERI>::value EvalCreator;
      return EvalCreator::eval(this,new IntParamsVoid);
    }

    /** Return a 2-body evaluator that computes two-electron integrals specific
        to linear R12 methods.  According to the convention in the
        literature, "g" stands for electron repulsion integral, "r" for the
        integral of r12 operator, and "t" for the commutator
        integrals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrR12.
        Implementation for this kind of TwoBodyInt is
        optional.

        @tparam NumCenters specifies the number of centers that carry basis functions.
        Valid values are 4, 3, and 2.
        */
    template <int NumCenters>
    DEPRECATED Ref< typename TwoBodyIntEvalType<NumCenters>::value > grt() {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::R12>::value EvalCreator;
      return EvalCreator::eval(this,new IntParamsVoid);
    }
    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use Gaussian geminals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrG12.
        Implementation for this kind of TwoBodyInt is optional. */
    template <int NumCenters>
    DEPRECATED Ref< typename TwoBodyIntEvalType<NumCenters>::value > g12(const Ref<IntParamsG12>& p) {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::G12>::value EvalCreator;
      return EvalCreator::eval(this,p);
    }
    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use Gaussian geminals.
        This particular implementation does not produce commutator integrals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrG12NC.
        Implementation for this kind of TwoBodyInt is optional. */
    template <int NumCenters>
    DEPRECATED Ref< typename TwoBodyIntEvalType<NumCenters>::value > g12nc(const Ref<IntParamsG12>& p) {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::G12NC>::value EvalCreator;
      return EvalCreator::eval(this,p);
    }
    /** Return a TwoBodyInt that computes two-electron integrals specific
        to relativistic explicitly correlated methods which use Gaussian geminals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrG12DKH.
        Implementation for this kind of TwoBodyInt is optional. */
    template <int NumCenters>
    Ref< typename TwoBodyIntEvalType<NumCenters>::value > g12dkh(const Ref<IntParamsG12>& p) {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::G12DKH>::value EvalCreator;
      return EvalCreator::eval(this,p);
    }

    /** Return the evaluator of two-body integrals with kernel
        \f$
          r_{12}^k g_{12}, \, k=-1,0,
        \f$
        where \f$ g_{12} \f$ is a geminal described by the IntParamsG12 object.
        These integrals are, for example, necessary in
        explicitly correlated methods which use Gaussian geminals.

        The evaluator will produce a set of integrals described by TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::R12_0_G12> for k=0
        and TwoBodyNCenterIntDescr<4,TwoBodyOperSet::R12_m1_G12> for k=-1.

        @tparam NumCenters specifies the number of centers that carry basis functions.
        Valid values are 4, 3, and 2.

        @note Implementation of this function is optional. The default implementation will throw FeatureNotImplemented .
      */
    template <int NumCenters>
    Ref< typename TwoBodyIntEvalType<NumCenters>::value > r12_k_g12(const Ref<IntParamsG12>& p, int k) {
      if (k == 0) {
        typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::R12_0_G12>::value EvalCreator;
        return EvalCreator::eval(this,p);
      }
      if (k == -1) {
        typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::R12_m1_G12>::value EvalCreator;
        return EvalCreator::eval(this,p);
      }
      std::ostringstream oss;
      oss << "Integral::r12_k_g12 can only be computed for k=-1 and 0 but k = " << k;
      throw FeatureNotImplemented(oss.str().c_str(),
                                  __FILE__, __LINE__);
    }

    /** Return the evaluator of two-body integrals with kernel
        \f$
          [g_{12},[\hat{T}_1,g_{12}]]
        \f$
        where \f$ g_{12} \f$ is a geminal described by the IntParamsG12 object.
        These integrals are, for example, necessary in
        explicitly correlated methods which use Gaussian geminals.

        The evaluator will produce a set of integrals described by TwoBodyNCenterIntDescr<NumCenters,TwoBodyOperSet::G12_T1_G12>.

        @tparam NumCenters specifies the number of centers that carry basis functions.
        Valid values are 4, 3, and 2.

        @note Implementation of this function is optional. The default implementation will throw FeatureNotImplemented .
      */
    template <int NumCenters>
    Ref< typename TwoBodyIntEvalType<NumCenters>::value > g12t1g12(const Ref<IntParamsG12>& p) {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::G12_T1_G12>::value EvalCreator;
      return EvalCreator::eval(this,p);
    }

    /** Return the evaluator of two-body integrals with kernel
        \f$
          \delta_3({\bf r}_1 - {\bf r}_2),
        \f$
        i.e. a one-electron overlap.

        The evaluator will produce a set of integrals described by TwoBodyIntDescrDelta.

        @tparam NumCenters specifies the number of centers that carry basis functions.
        Valid values are 4, 3, and 2.

        @note Implementation of this function is optional. The default implementation will throw sc::FeatureNotImplemented. Implemented in sc::IntegralLibint2
      */
    template <int NumCenters>
    Ref< typename TwoBodyIntEvalType<NumCenters>::value > delta_function() {
      typedef typename detail::EvalCreator<NumCenters,TwoBodyOperSet::DeltaFunction>::value EvalCreator;
      return EvalCreator::eval(this, new IntParamsVoid);
    }

  private:
    template <int NumCenters> friend struct sc::detail::ERIEvalCreator;
    template <int NumCenters> friend struct sc::detail::R12EvalCreator;
    template <int NumCenters> friend struct sc::detail::G12EvalCreator;
    template <int NumCenters> friend struct sc::detail::G12NCEvalCreator;
    template <int NumCenters> friend struct sc::detail::G12DKHEvalCreator;
    template <int NumCenters> friend struct sc::detail::R120G12EvalCreator;
    template <int NumCenters> friend struct sc::detail::R12m1G12EvalCreator;
    template <int NumCenters> friend struct sc::detail::G12T1G12EvalCreator;
    template <int NumCenters> friend struct sc::detail::DeltaFunctionEvalCreator;

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to linear R12 methods.  According to the convention in the
        literature, "g" stands for electron repulsion integral, "r" for the
        integral of r12 operator, and "t" for the commutator
        integrals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrR12.
        Implementation for this kind of TwoBodyInt is
        optional. */
    virtual Ref<TwoBodyInt> grt_4() INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> grt_3() INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> grt_2() INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use Gaussian geminals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrG12.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> g12_4(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> g12_3(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> g12_2(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to explicitly correlated methods which use Gaussian geminals.
	    This particular implementation does not produce commutator integrals.
	    This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrG12NC.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> g12nc_4(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> g12nc_3(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> g12nc_2(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals specific
        to relativistic explicitly correlated methods which use Gaussian geminals.
        This TwoBodyInt will produce a set of integrals described by TwoBodyIntDescrG12DKH.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> g12dkh_4(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> g12dkh_3(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> g12dkh_2(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals of TwoBodyOper::r12_0_g12.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> r120g12_4(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> r120g12_3(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> r120g12_2(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals of TwoBodyOper::r12_m1_g12.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> r12m1g12_4(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> r12m1g12_3(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> r12m1g12_2(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals of TwoBodyOper::g12t1g12.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> g12t1g12_4(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> g12t1g12_3(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> g12t1g12_2(const Ref<IntParamsG12>&) INTEGRALTHROWDEFAULT;

    /** Return a TwoBodyInt that computes two-electron integrals of TwoBodyOper::delta.
        Implementation for this kind of TwoBodyInt is optional. */
    virtual Ref<TwoBodyInt> delta_function_4() INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyThreeCenterInt> delta_function_3() INTEGRALTHROWDEFAULT;
    virtual Ref<TwoBodyTwoCenterInt> delta_function_2() INTEGRALTHROWDEFAULT;

};

// auxiliary type functions
namespace detail {

  template <> struct ERIEvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->electron_repulsion();
    }
  };
  template <> struct ERIEvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->electron_repulsion3();
    }
  };
  template <> struct ERIEvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->electron_repulsion2();
    }
  };

  template <> struct R12EvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->grt_4();
    }
  };
  template <> struct R12EvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->grt_3();
    }
  };
  template <> struct R12EvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->grt_2();
    }
  };

  template <> struct G12EvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12_4(params);
    }
  };
  template <> struct G12EvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12_3(params);
    }
  };
  template <> struct G12EvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12_2(params);
    }
  };

  template <> struct G12NCEvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12nc_4(params);
    }
  };
  template <> struct G12NCEvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12nc_3(params);
    }
  };
  template <> struct G12NCEvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12nc_2(params);
    }
  };

  template <> struct G12DKHEvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12dkh_4(params);
    }
  };
  template <> struct G12DKHEvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12dkh_3(params);
    }
  };
  template <> struct G12DKHEvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12dkh_2(params);
    }
  };

  template <> struct R120G12EvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->r120g12_4(params);
    }
  };
  template <> struct R120G12EvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->r120g12_3(params);
    }
  };
  template <> struct R120G12EvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->r120g12_2(params);
    }
  };

  template <> struct R12m1G12EvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->r12m1g12_4(params);
    }
  };
  template <> struct R12m1G12EvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->r12m1g12_3(params);
    }
  };
  template <> struct R12m1G12EvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->r12m1g12_2(params);
    }
  };

  template <> struct G12T1G12EvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12t1g12_4(params);
    }
  };
  template <> struct G12T1G12EvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12t1g12_3(params);
    }
  };
  template <> struct G12T1G12EvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsG12>& params) {
      return factory->g12t1g12_2(params);
    }
  };

  template <> struct DeltaFunctionEvalCreator<4> {
    static Ref< TwoBodyIntEvalType<4>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->delta_function_4();
    }
  };
  template <> struct DeltaFunctionEvalCreator<3> {
    static Ref< TwoBodyIntEvalType<3>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->delta_function_3();
    }
  };
  template <> struct DeltaFunctionEvalCreator<2> {
    static Ref< TwoBodyIntEvalType<2>::value >
    eval(Integral* factory, const Ref<IntParamsVoid>& params) {
      return factory->delta_function_2();
    }
  };

  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::ERI> {
    typedef ERIEvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::R12> {
    typedef R12EvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::G12> {
    typedef G12EvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::G12NC> {
    typedef G12NCEvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::G12DKH> {
    typedef G12DKHEvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::R12_0_G12> {
    typedef R120G12EvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::R12_m1_G12> {
    typedef R12m1G12EvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::G12_T1_G12> {
    typedef G12T1G12EvalCreator<NumCenters> value;
  };
  template <int NumCenters> struct EvalCreator<NumCenters,TwoBodyOperSet::DeltaFunction> {
    typedef DeltaFunctionEvalCreator<NumCenters> value;
  };

};

/// @}
// end of addtogroup ChemistryBasisIntegral

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
