//
// intv3.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <math/topology/pointbag.h>
#include <chemistry/qc/basis/integral.h>

class IntegralV3 : public Integral {
#   define CLASSNAME IntegralV3
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    IntegralV3(const RefGaussianBasisSet &b1=0,
               const RefGaussianBasisSet &b2=0,
               const RefGaussianBasisSet &b3=0,
               const RefGaussianBasisSet &b4=0);
    IntegralV3(StateIn&);
    IntegralV3(const RefKeyVal&);

    void save_data_state(StateOut&);
    
    CartesianIter * new_cartesian_iter(int);
    RedundantCartesianIter * new_redundant_cartesian_iter(int);
    RedundantCartesianSubIter * new_redundant_cartesian_sub_iter(int);
    SphericalTransformIter * new_spherical_transform_iter(int, int=0);
    
    RefOneBodyInt overlap();

    RefOneBodyInt kinetic();

    RefOneBodyInt point_charge(const RefPointChargeData& =0);

    RefOneBodyInt nuclear();

    RefOneBodyInt hcore();

    RefOneBodyInt efield_dot_vector(const RefEfieldDotVectorData& =0);

    RefOneBodyInt dipole(const RefDipoleData& =0);

    RefOneBodyDerivInt overlap_deriv();
                                     
    RefOneBodyDerivInt kinetic_deriv();
                                     
    RefOneBodyDerivInt nuclear_deriv();
                                     
    RefOneBodyDerivInt hcore_deriv();
                                     
    RefTwoBodyInt electron_repulsion();

    RefTwoBodyDerivInt electron_repulsion_deriv();
};


#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
