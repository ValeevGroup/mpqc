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
#include <util/misc/ccaenv.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/intcca/intcca.h>
#include <chemistry/qc/intcca/obintcca.h>
#include <chemistry/qc/intcca/tbintcca.h>

using namespace std;
using namespace sc;
using namespace Chemistry::QC::GaussianBasis;

static ClassDesc IntegralCCA_cd(
  typeid(IntegralCCA),"IntegralCCA",1,"public Integral",
  0, create<IntegralCCA>, create<IntegralCCA>);

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

IntegralCCA::IntegralCCA(const Ref<KeyVal> &keyval):
  Integral(keyval)
{

  initialize_transforms();

  string buffer = keyval->stringvalue("integral_buffer");
  if ( keyval->error() != KeyVal::OK ) buffer = "opaque";
  if ( buffer == "opaque" ) use_opaque_ = 1;
  else if ( buffer == "array" ) use_opaque_ = 0;
  else {
    ExEnv::err0() << indent << "unrecognized integral buffer type" << endl;
    abort();
  }

  factory_type_ = keyval->stringvalue("evaluator_factory");
  if ( keyval->error() != KeyVal::OK ) {
    factory_type_ = string("MPQC.IntegralEvaluatorFactory");
  }
  package_ = keyval->stringvalue("integral_package");
  if ( keyval->error() != KeyVal::OK ) {
    package_ = string("cints");
  }

  sc_molecule_ << keyval->describedclassvalue("molecule");
  if (sc_molecule_.null()) {
      ExEnv::err0() << indent << "molecule is required" << endl;
      abort();
  }

  gov::cca::Services &services = *CCAEnv::get_services();
  gov::cca::ports::BuilderService &bs = *CCAEnv::get_builder_service();
  gov::cca::TypeMap &type_map = *CCAEnv::get_type_map();
  gov::cca::ComponentID &my_id = *CCAEnv::get_component_id();

  // get eval factory
  fac_id_ = bs.createInstance("evaluator_factory",factory_type_,type_map);
  services.registerUsesPort("IntegralEvaluatorFactory",
                            "Chemistry.QC.GaussianBasis.IntegralEvaluatorFactory",
                            type_map);
  fac_con_ = bs.connect(my_id,"IntegralEvaluatorFactory",
                        fac_id_,"IntegralEvaluatorFactory");
  eval_factory_ = services.getPort("IntegralEvaluatorFactory");

  // set molecule on factory
  molecule_ = Chemistry::Chemistry_Molecule::_create();
  molecule_.initialize(sc_molecule_->natom(),"bohr");
  for( int i=0; i<sc_molecule_->natom(); ++i ) {
    molecule_.set_atomic_number( i, sc_molecule_->Z(i) );
    for( int j=0; j<3; ++j ) {
      molecule_.set_cart_coor( i, j, sc_molecule_->r(i,j) );
    }
  }
  eval_factory_.set_molecule(molecule_);

  // set package
  eval_factory_.set_integral_package(package_);

}

IntegralCCA::IntegralCCA(StateIn& s) :
  Integral(s)
{
  initialize_transforms();
}

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
  return new IntegralCCA(eval_factory_,use_opaque_);
}

CartesianIter *
IntegralCCA::new_cartesian_iter(int l)
{
  //return new CartesianIterV3(l);
  return new CartesianIterCints(l);
}

RedundantCartesianIter *
IntegralCCA::new_redundant_cartesian_iter(int l)
{
  //return new RedundantCartesianIterV3(l);
  return new RedundantCartesianIterCints(l);
}

RedundantCartesianSubIter *
IntegralCCA::new_redundant_cartesian_sub_iter(int l)
{
  //return new RedundantCartesianSubIterV3(l);
  return new RedundantCartesianSubIterCints(l);
}

SphericalTransformIter *
IntegralCCA::new_spherical_transform_iter(int l, int inv, int subl)
{

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
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, &Int1eCCA::overlap, 
                           use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::kinetic()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, &Int1eCCA::kinetic, 
                           use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::nuclear()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, &Int1eCCA::nuclear, 
                           use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::hcore()
{
  return new OneBodyIntCCA(this, bs1_, bs2_, eval_factory_, &Int1eCCA::hcore, 
                           use_opaque_ );
}

Ref<OneBodyInt>
IntegralCCA::point_charge(const Ref<PointChargeData>& dat)
{
//   return new PointChargeIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCCA::efield_dot_vector(const Ref<EfieldDotVectorData>&dat)
{
//   return new EfieldDotVectorIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCCA::dipole(const Ref<DipoleData>& dat)
{
//   return new DipoleIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCCA::quadrupole(const Ref<DipoleData>& dat)
{
//   throw std::runtime_error("IntegralV3 cannot compute quadrupole moment integrals yet. Try IntegralCints instead.");
}

Ref<OneBodyDerivInt>
IntegralCCA::overlap_deriv()
{
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::overlap_1der);
}

Ref<OneBodyDerivInt>
IntegralCCA::kinetic_deriv()
{
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::kinetic_1der);
}

Ref<OneBodyDerivInt>
IntegralCCA::nuclear_deriv()
{
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::nuclear_1der);
}

Ref<OneBodyDerivInt>
IntegralCCA::hcore_deriv()
{
//   return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::hcore_1der);
}

Ref<TwoBodyInt>
IntegralCCA::electron_repulsion()
{
  return new TwoBodyIntCCA(this, bs1_, bs2_, bs3_, bs4_, 
                           storage_, eval_factory_, use_opaque_ );
}

Ref<TwoBodyDerivInt>
IntegralCCA::electron_repulsion_deriv()
{
//   return new TwoBodyDerivIntV3(this, bs1_, bs2_, bs3_, bs4_, storage_);
}

void
IntegralCCA::set_basis(const Ref<GaussianBasisSet> &b1,
                       const Ref<GaussianBasisSet> &b2,
                       const Ref<GaussianBasisSet> &b3,
                       const Ref<GaussianBasisSet> &b4)
{
  free_transforms();
  Integral::set_basis(b1,b2,b3,b4);
  initialize_transforms();
}

void
IntegralCCA::free_transforms()
{

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
