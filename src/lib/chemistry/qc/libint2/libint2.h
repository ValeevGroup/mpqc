//
// libint2.h
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

// these provide integrals using the CINTS/libint routines

#ifndef _chemistry_qc_libint2_libint2_h
#define _chemistry_qc_libint2_libint2_h

#include <chemistry/qc/basis/integral.h>

namespace sc {

class SphericalTransformLibint2;
class ISphericalTransformLibint2;

/** IntegralLibint2 computes integrals between Gaussian basis functions. */
class IntegralLibint2 : public Integral {
  private:
    int maxl_;
    SphericalTransformLibint2 ***st_;
    ISphericalTransformLibint2 ***ist_;

    void free_transforms();
    void initialize_transforms();

    // Check if fully general contractions are present in any of the basis sets
    void check_fullgencon() const;

  public:
    IntegralLibint2(const Ref<GaussianBasisSet> &b1=0,
		  const Ref<GaussianBasisSet> &b2=0,
		  const Ref<GaussianBasisSet> &b3=0,
		  const Ref<GaussianBasisSet> &b4=0);
    IntegralLibint2(StateIn&);
    IntegralLibint2(const Ref<KeyVal>&);
    ~IntegralLibint2();

    void save_data_state(StateOut&);

    Integral* clone();

    /// implements Integral::cartesian_ordering()
    CartesianOrdering cartesian_ordering() const;

    size_t storage_required_eri(const Ref<GaussianBasisSet> &b1,
				const Ref<GaussianBasisSet> &b2 = 0,
				const Ref<GaussianBasisSet> &b3 = 0,
				const Ref<GaussianBasisSet> &b4 = 0);
    size_t storage_required_g12(const Ref<GaussianBasisSet> &b1,
				const Ref<GaussianBasisSet> &b2 = 0,
				const Ref<GaussianBasisSet> &b3 = 0,
				const Ref<GaussianBasisSet> &b4 = 0);
    size_t storage_required_g12nc(const Ref<GaussianBasisSet> &b1,
				  const Ref<GaussianBasisSet> &b2 = 0,
				  const Ref<GaussianBasisSet> &b3 = 0,
				  const Ref<GaussianBasisSet> &b4 = 0);
    size_t storage_required_g12dkh(const Ref<GaussianBasisSet> &b1,
                  const Ref<GaussianBasisSet> &b2 = 0,
                  const Ref<GaussianBasisSet> &b3 = 0,
                  const Ref<GaussianBasisSet> &b4 = 0);

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

    Ref<OneBodyInt> p4();

    Ref<OneBodyInt> hcore();

    Ref<OneBodyInt> efield_dot_vector(const Ref<EfieldDotVectorData>& =0);

    Ref<OneBodyInt> dipole(const Ref<DipoleData>& =0);

    Ref<OneBodyInt> quadrupole(const Ref<DipoleData>& =0);

    Ref<OneBodyDerivInt> overlap_deriv();

    Ref<OneBodyDerivInt> kinetic_deriv();

    Ref<OneBodyDerivInt> nuclear_deriv();

    Ref<OneBodyDerivInt> hcore_deriv();

    Ref<TwoBodyInt> electron_repulsion();
    Ref<TwoBodyThreeCenterInt> electron_repulsion3();
    Ref<TwoBodyTwoCenterInt> electron_repulsion2();

    Ref<TwoBodyDerivInt> electron_repulsion_deriv();

    Ref<TwoBodyInt> g12_4(const Ref<IntParamsG12>& p);

    Ref<TwoBodyInt> g12nc_4(const Ref<IntParamsG12>& p);
    Ref<TwoBodyThreeCenterInt> g12nc_3(const Ref<IntParamsG12>& p);
    Ref<TwoBodyTwoCenterInt> g12nc_2(const Ref<IntParamsG12>& p);

    Ref<TwoBodyInt> g12dkh_4(const Ref<IntParamsG12>& p);

    Ref<TwoBodyInt> r120g12_4(const Ref<IntParamsG12>& p);
    Ref<TwoBodyThreeCenterInt> r120g12_3(const Ref<IntParamsG12>& p);
    Ref<TwoBodyTwoCenterInt> r120g12_2(const Ref<IntParamsG12>& p);

    Ref<TwoBodyInt> r12m1g12_4(const Ref<IntParamsG12>& p);
    Ref<TwoBodyThreeCenterInt> r12m1g12_3(const Ref<IntParamsG12>& p);
    Ref<TwoBodyTwoCenterInt> r12m1g12_2(const Ref<IntParamsG12>& p);

    Ref<TwoBodyInt> g12t1g12_4(const Ref<IntParamsG12>& p);
    Ref<TwoBodyThreeCenterInt> g12t1g12_3(const Ref<IntParamsG12>& p);
    Ref<TwoBodyTwoCenterInt> g12t1g12_2(const Ref<IntParamsG12>& p);

    Ref<TwoBodyInt> delta_function_4();
    Ref<TwoBodyThreeCenterInt> delta_function_3();
    Ref<TwoBodyTwoCenterInt> delta_function_2();

    void set_basis(const Ref<GaussianBasisSet> &b1,
                   const Ref<GaussianBasisSet> &b2 = 0,
                   const Ref<GaussianBasisSet> &b3 = 0,
                   const Ref<GaussianBasisSet> &b4 = 0);
};

/* Libint2StaticInterface is an initializer class for the static part
   of libint's interface (one per executable) */
class Libint2StaticInterface {
    bool ready;

    public:
    Libint2StaticInterface();
    ~Libint2StaticInterface() { ready = false; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
