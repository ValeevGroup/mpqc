//
// intcca.h
//
// Copyright (C) 2004 Sandia National Laboratories
//
// Author: Joe Kenny <jpkenny@sandia.gov>
// Maintainer: Joe Kenny
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

// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_intcca_intcca_h
#define _chemistry_qc_intcca_intcca_h

#include <gov_cca.hh>
#include <chemistry/qc/basis/integral.h>
#include <Chemistry_QC_GaussianBasis_IntegralEvaluatorFactory.hh>
#include <Chemistry_Chemistry_Molecule.hh>
#include <chemistry/molecule/molecule.h>

// INTV3 includes
//#include <chemistry/qc/intv3/cartitv3.h>
//#include <chemistry/qc/intv3/tformv3.h>

// CINTS includes
#include <chemistry/qc/cints/cartit.h>
#include <chemistry/qc/cints/tform.h>

using namespace Chemistry::QC::GaussianBasis;

namespace sc {

/** IntegralCCA provides an MPQC client for CCA IntegralEvaluator components. */
class IntegralCCA : public Integral {
  private:
    int maxl_;
    bool use_opaque_;
    gov::cca::ComponentID fac_id_;
    gov::cca::ConnectionID fac_con_;
    Ref<Molecule> sc_molecule_;
    Chemistry::Chemistry_Molecule molecule_;
    std::string factory_type_;
    std::string package_;

    // INTV3 verion
//    SphericalTransformV3 ***st_;
//    ISphericalTransformV3 ***ist_;

    // CINTS version
    SphericalTransformCints ***st_;
    ISphericalTransformCints ***ist_;

    void free_transforms();
    void initialize_transforms();
    IntegralEvaluatorFactory eval_factory_;

  public:

    /** This constructor is used when the framework is not embedded. */
    IntegralCCA(IntegralEvaluatorFactory eval_factory, bool use_opaque,
                const Ref<GaussianBasisSet> &b1=0,
                const Ref<GaussianBasisSet> &b2=0,
                const Ref<GaussianBasisSet> &b3=0,
                const Ref<GaussianBasisSet> &b4=0);

    IntegralCCA(StateIn&);

    /** The KeyVal constructor.
        This constructor is used when the framework is embedded.
        The following keywords are read:

        <dl>
        <dt><tt>evaluator_factory</tt><dd> This gives the symbol name of a 
        CCA IntegralEvaluatorFactory component.  This symbol name should
        also appear in the cca-load argument.  The default is
        <tt>MPQC.IntegralEvaluatorFactory</tt>.

        <dt><tt>integral_package</tt><dd> This gives the name of the
        integrals package to use (<tt>intv3</tt> or <tt>cints</tt>).
        The default is <tt>cints</tt>.

        <dt><tt>molecule</tt><dd> This gives a molecule object, it is required.
        </dl>
    */

    IntegralCCA(const Ref<KeyVal>&);

    ~IntegralCCA();

    void save_data_state(StateOut&);

    Integral* clone();
   
    CartesianIter * new_cartesian_iter(int);
    RedundantCartesianIter * new_redundant_cartesian_iter(int);
    RedundantCartesianSubIter * new_redundant_cartesian_sub_iter(int);
    SphericalTransformIter * new_spherical_transform_iter(int l,
                                                          int inv=0,
                                                          int subl=-1);
    const SphericalTransform * spherical_transform(int l,
                                                   int inv=0, int subl=-1);

    Ref<OneBodyInt> overlap();

    Ref<OneBodyInt> kinetic();

    Ref<OneBodyInt> point_charge(const Ref<PointChargeData>& =0);

    Ref<OneBodyInt> nuclear();

    Ref<OneBodyInt> hcore();

    Ref<OneBodyInt> efield_dot_vector(const Ref<EfieldDotVectorData>& =0);

    Ref<OneBodyInt> dipole(const Ref<DipoleData>& =0);

    Ref<OneBodyInt> quadrupole(const Ref<DipoleData>& =0);

    Ref<OneBodyDerivInt> overlap_deriv();
                                     
    Ref<OneBodyDerivInt> kinetic_deriv();
                                     
    Ref<OneBodyDerivInt> nuclear_deriv();
                                     
    Ref<OneBodyDerivInt> hcore_deriv();
                                     
    Ref<TwoBodyInt> electron_repulsion();

    Ref<TwoBodyDerivInt> electron_repulsion_deriv();

    void set_basis(const Ref<GaussianBasisSet> &b1,
                   const Ref<GaussianBasisSet> &b2 = 0,
                   const Ref<GaussianBasisSet> &b3 = 0,
                   const Ref<GaussianBasisSet> &b4 = 0);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
