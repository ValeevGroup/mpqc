//
// cints.cc
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

#include <stdexcept>

#include <util/state/stateio.h>
#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <chemistry/qc/cints/cints.h>
#include <chemistry/qc/cints/cartit.h>
#include <chemistry/qc/cints/tform.h>
#include <chemistry/qc/cints/obintcints.h>
#include <chemistry/qc/cints/tbintcints.h>
#include <chemistry/qc/cints/eri.h>
#include <chemistry/qc/cints/grt.h>

using namespace std;
using namespace sc;

inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

static ClassDesc IntegralCints_cd(
  typeid(IntegralCints),"IntegralCints",1,"public Integral",
  0, create<IntegralCints>, create<IntegralCints>);

IntegralCints::IntegralCints(const Ref<GaussianBasisSet> &b1,
			     const Ref<GaussianBasisSet> &b2,
			     const Ref<GaussianBasisSet> &b3,
			     const Ref<GaussianBasisSet> &b4):
  Integral(b1,b2,b3,b4)
{
  initialize_transforms();
}

IntegralCints::IntegralCints(StateIn& s) :
  Integral(s)
{
  initialize_transforms();
}

IntegralCints::IntegralCints(const Ref<KeyVal>& k) :
  Integral(k)
{
  initialize_transforms();
}

void
IntegralCints::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

IntegralCints::~IntegralCints()
{
  free_transforms();
}

Integral*
IntegralCints::clone()
{
  return new IntegralCints;
}

size_t
IntegralCints::storage_required_eri(const Ref<GaussianBasisSet> &b1,
				    const Ref<GaussianBasisSet> &b2,
				    const Ref<GaussianBasisSet> &b3,
				    const Ref<GaussianBasisSet> &b4)
{
  return EriCints::storage_required(b1,b2,b3,b4);
}

size_t
IntegralCints::storage_required_grt(const Ref<GaussianBasisSet> &b1,
				    const Ref<GaussianBasisSet> &b2,
				    const Ref<GaussianBasisSet> &b3,
				    const Ref<GaussianBasisSet> &b4)
{
  return GRTCints::storage_required(b1,b2,b3,b4);
}

CartesianIter *
IntegralCints::new_cartesian_iter(int l)
{
  return new CartesianIterCints(l);
}

RedundantCartesianIter *
IntegralCints::new_redundant_cartesian_iter(int l)
{
  return new RedundantCartesianIterCints(l);
}

RedundantCartesianSubIter *
IntegralCints::new_redundant_cartesian_sub_iter(int l)
{
  return new RedundantCartesianSubIterCints(l);
}

SphericalTransformIter *
IntegralCints::new_spherical_transform_iter(int l, int inv, int subl)
{
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
IntegralCints::spherical_transform(int l, int inv, int subl)
{
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
IntegralCints::overlap()
{
  return new OneBodyIntCints(this, bs1_, bs2_, &Int1eCints::overlap);
}

Ref<OneBodyInt>
IntegralCints::kinetic()
{
  return new OneBodyIntCints(this, bs1_, bs2_, &Int1eCints::kinetic);
}

Ref<OneBodyInt>
IntegralCints::nuclear()
{
  return new OneBodyIntCints(this, bs1_, bs2_, &Int1eCints::nuclear);
}

Ref<OneBodyInt>
IntegralCints::p4()
{
  throw FeatureNotImplemented("IntegralCints cannot compute p4 integrals. Try other integral factories.");
}

Ref<OneBodyInt>
IntegralCints::hcore()
{
  return new OneBodyIntCints(this, bs1_, bs2_, &Int1eCints::hcore);
}

Ref<OneBodyInt>
IntegralCints::point_charge(const Ref<PointChargeData>& dat)
{
  ExEnv::errn() << scprintf("IntegralCints::point_charge() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new PointChargeIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCints::efield_dot_vector(const Ref<EfieldDotVectorData>&dat)
{
  ExEnv::errn() << scprintf("IntegralCints::efield_dot_vector() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new EfieldDotVectorIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralCints::dipole(const Ref<DipoleData>& dat)
{
  Ref<OneBodyIntCints> dipoleint = new OneBodyIntCints(this, bs1_, bs2_, &Int1eCints::edipole);
  dipoleint->set_multipole_origin(dat);
  return dipoleint;
}

Ref<OneBodyInt>
IntegralCints::quadrupole(const Ref<DipoleData>& dat)
{
  Ref<OneBodyIntCints> quadint = new OneBodyIntCints(this, bs1_, bs2_, &Int1eCints::equadrupole);
  quadint->set_multipole_origin(dat);
  return quadint;
}

Ref<OneBodyDerivInt>
IntegralCints::overlap_deriv()
{
  ExEnv::errn() << scprintf("IntegralCints::overlap_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::overlap_1der);
}

Ref<OneBodyDerivInt>
IntegralCints::kinetic_deriv()
{
  ExEnv::errn() << scprintf("IntegralCints::kinetic_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::kinetic_1der);
}

Ref<OneBodyDerivInt>
IntegralCints::nuclear_deriv()
{
  ExEnv::errn() << scprintf("IntegralCints::nuclear_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::nuclear_1der);
}

Ref<OneBodyDerivInt>
IntegralCints::hcore_deriv()
{
  ExEnv::errn() << scprintf("IntegralCints::hcore_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::hcore_1der);
}

Ref<TwoBodyInt>
IntegralCints::electron_repulsion()
{
  return new TwoBodyIntCints(this, bs1_, bs2_, bs3_, bs4_, storage_, TwoBodyOperSet::ERI);
}

Ref<TwoBodyInt>
IntegralCints::grt_4()
{
  return new TwoBodyIntCints(this, bs1_, bs2_, bs3_, bs4_, storage_, TwoBodyOperSet::R12);
}

namespace {
  int maxl(const Ref<GaussianBasisSet>& bs) {
    int result = -1;
    if (bs.nonnull())
      result = bs->max_angular_momentum();
    return result;
  }
}

void
IntegralCints::set_basis(const Ref<GaussianBasisSet> &b1,
			 const Ref<GaussianBasisSet> &b2,
			 const Ref<GaussianBasisSet> &b3,
			 const Ref<GaussianBasisSet> &b4)
{
  Integral::set_basis(b1,b2,b3,b4);
  check_fullgencon();
  const int maxl_new = std::max(std::max(maxl(b1),maxl(b2)),
                                std::max(maxl(b3),maxl(b4)) );
  if (maxl_new > maxl_) {
    free_transforms();
    initialize_transforms();
  }
}

void
IntegralCints::free_transforms()
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
  if (maxl_ >= 0) {
    delete[] st_;
    delete[] ist_;
  }
  st_ = NULL;
  ist_ = NULL;
}

void
IntegralCints::initialize_transforms()
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


static bool has_fullgencon(const Ref<GaussianBasisSet>&);

void
IntegralCints::check_fullgencon() const
{
  if ( has_fullgencon(bs1_) ||
       has_fullgencon(bs2_) ||
       has_fullgencon(bs3_) ||
       has_fullgencon(bs4_) ) {
      throw std::runtime_error("IntegralCints cannot handle basis sets with fully general contractions yet, try IntegralV3 instead");
  }
}

bool
has_fullgencon(const Ref<GaussianBasisSet>& bs)
{
  bool has_it = false;
  int nshell = bs->nshell();
  for(int i=0; i<nshell; i++) {
    GaussianShell& shell = bs->shell(i);
    int minam = shell.min_am();
    int maxam = shell.max_am();
    if (minam != maxam)
      has_it = true;
  }

  return has_it;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
