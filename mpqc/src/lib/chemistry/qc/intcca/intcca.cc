//
// intcca.cc
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

#include <stdexcept>

#include <util/state/stateio.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/intcca/intcca.h>
#include <chemistry/qc/intcca/obintcca.h>
#include <chemistry/qc/intcca/tbintcca.h>

using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;

static ClassDesc IntegralCCA_cd(
  typeid(IntegralCCA),"IntegralCCA",1,"public Integral",
  0, 0, create<IntegralCCA>);

extern Ref<Integral> default_integral;

/* may need to add optional "eval_factory" argument to this method in
integral class to get this capability
Integral*
Integral::get_default_integral()
{
  if (default_integral.null())
    default_integral = new IntegralCCA();

  return default_integral;
}
*/

IntegralCCA::IntegralCCA(IntegralEvaluatorFactory eval_factory, bool use_opaque,
                         const Ref<GaussianBasisSet> &b1,
                         const Ref<GaussianBasisSet> &b2,
                         const Ref<GaussianBasisSet> &b3,
                         const Ref<GaussianBasisSet> &b4):
  Integral(b1,b2,b3,b4), eval_factory_(eval_factory)
{
  use_opaque_ = use_opaque;
  initialize_transforms();
}

IntegralCCA::IntegralCCA(StateIn& s) :
  Integral(s)
{
  initialize_transforms();
}

//IntegralCCA::IntegralCCA(const Ref<KeyVal>& k) :
//  Integral(k)
//{
//  initialize_transforms();
//}

void
IntegralCCA::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

IntegralCCA::~IntegralCCA()
{
  free_transforms();
}

Integral*
IntegralCCA::clone()
{
  //ExEnv::out0() << "IntegralCCA::clone() called\n";
  return new IntegralCCA(eval_factory_,use_opaque_);
}

CartesianIter *
IntegralCCA::new_cartesian_iter(int l)
{
  //ExEnv::out0() << "IntegralCCA::new_cartesian_iter() called\n";
  //return new CartesianIterV3(l);
  return new CartesianIterCints(l);
}

RedundantCartesianIter *
IntegralCCA::new_redundant_cartesian_iter(int l)
{
  //ExEnv::out0() << "IntegralCCA::new_redundant_cartesian_iter() called\n";
  //return new RedundantCartesianIterV3(l);
  return new RedundantCartesianIterCints(l);
}

RedundantCartesianSubIter *
IntegralCCA::new_redundant_cartesian_sub_iter(int l)
{
  //ExEnv::out0() << "IntegralCCA::new_redundant_cartesian_sub_iter() called\n";
  //return new RedundantCartesianSubIterV3(l);
  return new RedundantCartesianSubIterCints(l);
}

SphericalTransformIter *
IntegralCCA::new_spherical_transform_iter(int l, int inv, int subl)
{
  //ExEnv::out0() << "IntegralCCA::new_spherical_transform iter() called\n";

  // INTV3 version
/*  if (l>maxl_ || l<0) {
      ExEnv::errn() << "IntegralV3::new_spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      ExEnv::errn() << "IntegralV3::new_spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
    }
  return new SphericalTransformIter(st_[l][(l-subl)/2]);
*/
 
  // CINTS version
  if (l>maxl_ || l<0) {
      ExEnv::errn() << "IntegralCints::new_spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      ExEnv::errn() << "IntegralCints::new_spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
    }
  return new SphericalTransformIter(st_[l][(l-subl)/2]);


}

const SphericalTransform *
IntegralCCA::spherical_transform(int l, int inv, int subl)
{
  //ExEnv::out0() << "IntegralCCA::spherical_transform() called\n";	

  // INTV3 version
/*  if (l>maxl_ || l<0) {
      ExEnv::errn() << "IntegralV3::spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      ExEnv::errn() << "IntegralV3::spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return ist_[l][(l-subl)/2];
    }
  return st_[l][(l-subl)/2];
*/

  // CINTS version
  if (l>maxl_ || l<0) {
      ExEnv::errn() << "IntegralCints::spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      ExEnv::errn() << "IntegralCints::spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return ist_[l][(l-subl)/2];
    }
  return st_[l][(l-subl)/2];

}

Ref<OneBodyInt>
IntegralCCA::overlap()
{
  //ExEnv::out0() << "IntegralCCA::overlap() called\n";
  return new OneBodyIntCCA(this, bs1_, bs2_, &Int1eCCA::overlap, 
                           eval_factory_, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::kinetic()
{
  //ExEnv::out0() << "IntegralCCA::kinetic() called\n";
  return new OneBodyIntCCA(this, bs1_, bs2_, &Int1eCCA::kinetic, 
                           eval_factory_, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::nuclear()
{
  //ExEnv::out0() << "IntegralCCA::nuclear() called\n";
  return new OneBodyIntCCA(this, bs1_, bs2_, &Int1eCCA::nuclear, 
                           eval_factory_, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::hcore()
{
  //ExEnv::out0() << "IntegralCCA::hcore() called\n";
  return new OneBodyIntCCA(this, bs1_, bs2_, &Int1eCCA::hcore, 
                           eval_factory_, use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::point_charge(const Ref<PointChargeData>& dat)
{
  //ExEnv::out0() << "IntegralCCA::point_charge() (not implemented) called\n";
//   return new PointChargeIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCCA::efield_dot_vector(const Ref<EfieldDotVectorData>&dat)
{
  //ExEnv::out0() << "IntegralCCA::efield_dot_vector() (not implemented) called\n";
//   return new EfieldDotVectorIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCCA::dipole(const Ref<DipoleData>& dat)
{
//ExEnv::out0() << "IntegralCCA::dipole() (not implemented) called\n";
//   return new DipoleIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCCA::quadrupole(const Ref<DipoleData>& dat)
{
//ExEnv::out0() << "IntegralCCA::quadrupole() (not implemented) called\n";
//   throw std::runtime_error("IntegralV3 cannot compute quadrupole moment integrals yet. Try IntegralCints instead.");
}

Ref<OneBodyDerivInt>
IntegralCCA::overlap_deriv()
{
//ExEnv::out0() << "IntegralCCA::overlap_deriv() (not implemented) called\n";
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::overlap_1der);
}

Ref<OneBodyDerivInt>
IntegralCCA::kinetic_deriv()
{
//ExEnv::out0() << "IntegralCCA::kinetic_deriv() (not implemented) called\n";
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::kinetic_1der);
}

Ref<OneBodyDerivInt>
IntegralCCA::nuclear_deriv()
{
//ExEnv::out0() << "IntegralCCA::nuclear_deriv() (not implemented) called\n";
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::nuclear_1der);
}

Ref<OneBodyDerivInt>
IntegralCCA::hcore_deriv()
{
//ExEnv::out0() << "IntegralCCA::hcore_deriv() (not implemented) called\n";
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::hcore_1der);
}

Ref<TwoBodyInt>
IntegralCCA::electron_repulsion()
{
  //ExEnv::out0() << "IntegralCCA::electron_repulsion() called\n";
  return new TwoBodyIntCCA(this, bs1_, bs2_, bs3_, bs4_, 
                           storage_, eval_factory_, use_opaque_ );
}

Ref<TwoBodyDerivInt>
IntegralCCA::electron_repulsion_deriv()
{
//ExEnv::out0() << "IntegralCCA::electron_repulsion_deriv() (not implemented) called\n";
//   return new TwoBodyDerivIntV3(this, bs1_, bs2_, bs3_, bs4_, storage_);
}

void
IntegralCCA::set_basis(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2,
                       const Ref<GaussianBasisSet> &b3,
                       const Ref<GaussianBasisSet> &b4)
{
  //ExEnv::out0() << "IntegralCCA::set_basis() called\n";
  free_transforms();
  Integral::set_basis(b1,b2,b3,b4);
  initialize_transforms();
}

void
IntegralCCA::free_transforms()
{
  //ExEnv::out0() << "IntegralCCA::free_transforms() called\n";

  // INTV3 version
/*  int i,j;
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
*/

  // CINTS version
  int i,j;
  for (i=0; i<=maxl_; i++) {
      for (j=0; j<=i/2; j++) {
          delete st_[i][j];
          delete ist_[i][j];
        }
      delete[] st_[i];
      delete[] ist_[i];
    }
  if (maxl_ >= 0) {
    delete[] st_;
    delete[] ist_;
  }
  st_ = NULL;
  ist_ = NULL;

}

 void
 IntegralCCA::initialize_transforms()
 {
   //ExEnv::out0() << "IntegralCCA::initialize_transforms() called\n";
 
 
   // INTV3 version
/*   maxl_ = -1;
   int maxam;
   maxam = bs1_.nonnull()?bs1_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
   maxam = bs2_.nonnull()?bs2_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
   maxam = bs3_.nonnull()?bs3_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;
   maxam = bs4_.nonnull()?bs4_->max_angular_momentum():-1;
   if (maxl_ < maxam) maxl_ = maxam;

   st_ = new SphericalTransformV3**[maxl_+1];
   ist_ = new ISphericalTransformV3**[maxl_+1];;
   int i,j;
   for (i=0; i<=maxl_; i++) {
       st_[i] = new SphericalTransformV3*[i/2+1];
       ist_[i] = new ISphericalTransformV3*[i/2+1];
       for (j=0; j<=i/2; j++) {
           st_[i][j] = new SphericalTransformV3(i,i-2*j);
           ist_[i][j] = new ISphericalTransformV3(i,i-2*j);
         }
     }
*/


  // CINTS version
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

  if (maxl_ >= 0) {
    st_ = new SphericalTransformCints**[maxl_+1];
    ist_ = new ISphericalTransformCints**[maxl_+1];;
    int i,j;
    for (i=0; i<=maxl_; i++) {
      st_[i] = new SphericalTransformCints*[i/2+1];
      ist_[i] = new ISphericalTransformCints*[i/2+1];
      for (j=0; j<=i/2; j++) {
        st_[i][j] = new SphericalTransformCints(i,i-2*j);
        ist_[i][j] = new ISphericalTransformCints(i,i-2*j);
      }
    }
  }
  else {
    st_ = NULL;
    ist_ = NULL;
  }

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
