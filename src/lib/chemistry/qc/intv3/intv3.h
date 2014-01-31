//
// intv3.h
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

// these provide integrals using the libintv2 routines

#ifndef _chemistry_qc_intv3_intv3_h
#define _chemistry_qc_intv3_intv3_h

#include <chemistry/qc/basis/integral.h>

namespace sc {

class SphericalTransformV3;
class ISphericalTransformV3;

/** IntegralV3 computes integrals between Gaussian basis functions. */
class IntegralV3 : public Integral {
  private:
    int maxl_;
    SphericalTransformV3 ***st_;
    ISphericalTransformV3 ***ist_;

    void free_transforms();
    void initialize_transforms();
  public:
    IntegralV3(const Ref<GaussianBasisSet> &b1=0,
               const Ref<GaussianBasisSet> &b2=0,
               const Ref<GaussianBasisSet> &b3=0,
               const Ref<GaussianBasisSet> &b4=0);
    IntegralV3(StateIn&);
    IntegralV3(const Ref<KeyVal>&);
    ~IntegralV3();

    void save_data_state(StateOut&);

    Integral* clone();

    /// implements Integral::cartesian_ordering()
    CartesianOrdering cartesian_ordering() const { return IntV3CartesianOrdering; }
    
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

    Ref<OneBodyOneCenterInt> point_charge1(const Ref<PointChargeData>&);

    Ref<OneBodyInt> nuclear();

    Ref<OneBodyInt> p_dot_nuclear_p();

    Ref<OneBodyInt> p4();
    
    Ref<OneBodyInt> hcore();

    Ref<OneBodyInt> efield(const Ref<IntParamsOrigin>& =0);

    Ref<OneBodyInt> efield_dot_vector(const Ref<EfieldDotVectorData>& =0);

    Ref<OneBodyInt> dipole(const Ref<IntParamsOrigin>& =0);

    Ref<OneBodyInt> quadrupole(const Ref<IntParamsOrigin>& =0);

    Ref<OneBodyDerivInt> overlap_deriv();
                                     
    Ref<OneBodyDerivInt> kinetic_deriv();
                                     
    Ref<OneBodyDerivInt> nuclear_deriv();
                                     
    Ref<OneBodyDerivInt> hcore_deriv();
                                     
    Ref<TwoBodyInt> electron_repulsion();

    Ref<TwoBodyTwoCenterInt> electron_repulsion2();

    Ref<TwoBodyThreeCenterInt> electron_repulsion3();

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
