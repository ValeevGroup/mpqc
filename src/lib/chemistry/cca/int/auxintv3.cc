//
// auxintv3.cc
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

#include <util/class/scexception.h>
#include <chemistry/cca/int/intcca.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/intv3/tformv3.h>

using namespace std;
using namespace sc;

CartesianIter *
IntegralCCA::new_cartesian_iterV3(int l)
{ 
  return new CartesianIterV3(l);
}

RedundantCartesianIter *
IntegralCCA::new_redundant_cartesian_iterV3(int l)
{
  return new RedundantCartesianIterV3(l);
}

RedundantCartesianSubIter *
IntegralCCA::new_redundant_cartesian_sub_iterV3(int l)
{
  return new RedundantCartesianSubIterV3(l);
}

SphericalTransformIter *
IntegralCCA::new_spherical_transform_iterV3(int l, int inv, int subl)
{
  if (l>maxl_ || l<0)
      throw ProgrammingError("new_spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("new_spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
  return new SphericalTransformIter(st_[l][(l-subl)/2]);
}

const SphericalTransform *
IntegralCCA::spherical_transformV3(int l, int inv, int subl)
{
  if (l>maxl_ || l<0)
      throw ProgrammingError("spherical_transform_iter: bad l",
                             __FILE__,__LINE__);
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0)
      throw ProgrammingError("spherical_transform_iter: bad subl",
                             __FILE__,__LINE__);
  if (inv)
      return ist_[l][(l-subl)/2];
  return st_[l][(l-subl)/2];
}

void
IntegralCCA::free_transformsV3()
{
  int i,j;
  for (i=0; i<=maxl_; i++) {
    for (j=0; j<=i/2; j++) {
        delete st_[i][j];
        delete ist_[i][j];
    }
    delete[] st_[i];
    delete[] ist_[i];
  }
  delete[] st_;
  delete[] ist_;
  st_ = 0;
  ist_ = 0;
}

void
IntegralCCA::initialize_transformsV3()
{
  maxl_ = -1;
  int maxam;
  maxam = bs1_.nonnull()?bs1_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;
  maxam = bs2_.nonnull()?bs2_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;
  maxam = bs3_.nonnull()?bs3_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;
  maxam = bs4_.nonnull()?bs4_->max_angular_momentum():-1;
  if (maxl_ < maxam) maxl_ = maxam;

  st_ = new SphericalTransform**[maxl_+1];
  ist_ = new ISphericalTransform**[maxl_+1];
  int i,j;
  for (i=0; i<=maxl_; i++) {
    st_[i] = new SphericalTransform*[i/2+1];
    ist_[i] = new ISphericalTransform*[i/2+1];
    for (j=0; j<=i/2; j++) {
      st_[i][j] = new SphericalTransformV3(i,i-2*j);
      ist_[i][j] = new ISphericalTransformV3(i,i-2*j);
    }
  }

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
