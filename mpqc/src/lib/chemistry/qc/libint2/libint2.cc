//
// libint2.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/libint2/libint2.h>
#include <chemistry/qc/libint2/cartit.h>
#include <chemistry/qc/libint2/tform.h>
#include <chemistry/qc/libint2/obintlibint2.h>
#include <chemistry/qc/libint2/tbintlibint2.h>
#include <chemistry/qc/libint2/eri.h>
#include <chemistry/qc/libint2/g12.h>

using namespace std;
using namespace sc;

inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

static ClassDesc IntegralLibint2_cd(
  typeid(IntegralLibint2),"IntegralLibint2",1,"public Integral",
  0, create<IntegralLibint2>, create<IntegralLibint2>);

IntegralLibint2::IntegralLibint2(const Ref<GaussianBasisSet> &b1,
			     const Ref<GaussianBasisSet> &b2,
			     const Ref<GaussianBasisSet> &b3,
			     const Ref<GaussianBasisSet> &b4):
  Integral(b1,b2,b3,b4)
{
  initialize_transforms();
}

IntegralLibint2::IntegralLibint2(StateIn& s) :
  Integral(s)
{
  initialize_transforms();
}

IntegralLibint2::IntegralLibint2(const Ref<KeyVal>& k) :
  Integral(k)
{
  initialize_transforms();
}

void
IntegralLibint2::save_data_state(StateOut& s)
{
  Integral::save_data_state(s);
}

IntegralLibint2::~IntegralLibint2()
{
  free_transforms();
}

Integral*
IntegralLibint2::clone()
{
  return new IntegralLibint2;
}

size_t
IntegralLibint2::storage_required_eri(const Ref<GaussianBasisSet> &b1,
				    const Ref<GaussianBasisSet> &b2,
				    const Ref<GaussianBasisSet> &b3,
				    const Ref<GaussianBasisSet> &b4)
{
  return EriLibint2::storage_required(b1,b2,b3,b4);
}

size_t
IntegralLibint2::storage_required_g12(const Ref<GaussianBasisSet> &b1,
				     const Ref<GaussianBasisSet> &b2,
				     const Ref<GaussianBasisSet> &b3,
				     const Ref<GaussianBasisSet> &b4)
{
  return G12Libint2::storage_required(b1,b2,b3,b4);
}

size_t
IntegralLibint2::storage_required_grt(const Ref<GaussianBasisSet> &b1,
				    const Ref<GaussianBasisSet> &b2,
				    const Ref<GaussianBasisSet> &b3,
				    const Ref<GaussianBasisSet> &b4)
{
  throw FeatureNotImplemented("IntegralLibint2::storage_required_grt() -- IntegralLibint2 does not implement GRT integrals yet.", __FILE__, __LINE__);
}

CartesianIter *
IntegralLibint2::new_cartesian_iter(int l)
{
  return new CartesianIterLibint2(l);
}

RedundantCartesianIter *
IntegralLibint2::new_redundant_cartesian_iter(int l)
{
  return new RedundantCartesianIterLibint2(l);
}

RedundantCartesianSubIter *
IntegralLibint2::new_redundant_cartesian_sub_iter(int l)
{
  return new RedundantCartesianSubIterLibint2(l);
}

SphericalTransformIter *
IntegralLibint2::new_spherical_transform_iter(int l, int inv, int subl)
{
  if (l>maxl_ || l<0) {
      ExEnv::errn() << "IntegralLibint2::new_spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      ExEnv::errn() << "IntegralLibint2::new_spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return new SphericalTransformIter(ist_[l][(l-subl)/2]);
    }
  return new SphericalTransformIter(st_[l][(l-subl)/2]);
}

const SphericalTransform *
IntegralLibint2::spherical_transform(int l, int inv, int subl)
{
  if (l>maxl_ || l<0) {
      ExEnv::errn() << "IntegralLibint2::spherical_transform_iter: bad l" << endl;
      abort();
    }
  if (subl == -1) subl = l;
  if (subl < 0 || subl > l || (l-subl)%2 != 0) {
      ExEnv::errn() << "IntegralLibint2::spherical_transform_iter: bad subl" << endl;
      abort();
    }
  if (inv) {
      return ist_[l][(l-subl)/2];
    }
  return st_[l][(l-subl)/2];
}

Ref<OneBodyInt>
IntegralLibint2::overlap()
{
  return new OneBodyIntLibint2(this, bs1_, bs2_, &Int1eLibint2::overlap);
}

Ref<OneBodyInt>
IntegralLibint2::kinetic()
{
  return new OneBodyIntLibint2(this, bs1_, bs2_, &Int1eLibint2::kinetic);
}

Ref<OneBodyInt>
IntegralLibint2::nuclear()
{
  return new OneBodyIntLibint2(this, bs1_, bs2_, &Int1eLibint2::nuclear);
}

Ref<OneBodyInt>
IntegralLibint2::hcore()
{
  return new OneBodyIntLibint2(this, bs1_, bs2_, &Int1eLibint2::hcore);
}

Ref<OneBodyInt>
IntegralLibint2::point_charge(const Ref<PointChargeData>& dat)
{
  ExEnv::errn() << scprintf("IntegralLibint2::point_charge() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new PointChargeIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralLibint2::efield_dot_vector(const Ref<EfieldDotVectorData>&dat)
{
  ExEnv::errn() << scprintf("IntegralLibint2::efield_dot_vector() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new EfieldDotVectorIntV3(this, bs1_, bs2_, dat);
}

Ref<OneBodyInt>
IntegralLibint2::dipole(const Ref<DipoleData>& dat)
{
  Ref<OneBodyIntLibint2> dipoleint = new OneBodyIntLibint2(this, bs1_, bs2_, &Int1eLibint2::edipole);
  dipoleint->set_multipole_origin(dat);
  return dipoleint;
}

Ref<OneBodyInt>
IntegralLibint2::quadrupole(const Ref<DipoleData>& dat)
{
  Ref<OneBodyIntLibint2> quadint = new OneBodyIntLibint2(this, bs1_, bs2_, &Int1eLibint2::equadrupole);
  quadint->set_multipole_origin(dat);
  return quadint;
}

Ref<OneBodyDerivInt>
IntegralLibint2::overlap_deriv()
{
  ExEnv::errn() << scprintf("IntegralLibint2::overlap_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::overlap_1der);
}

Ref<OneBodyDerivInt>
IntegralLibint2::kinetic_deriv()
{
  ExEnv::errn() << scprintf("IntegralLibint2::kinetic_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::kinetic_1der);
}

Ref<OneBodyDerivInt>
IntegralLibint2::nuclear_deriv()
{
  ExEnv::errn() << scprintf("IntegralLibint2::nuclear_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::nuclear_1der);
}

Ref<OneBodyDerivInt>
IntegralLibint2::hcore_deriv()
{
  ExEnv::errn() << scprintf("IntegralLibint2::hcore_deriv() is not yet implemented.\n");
  ExEnv::errn() << scprintf("Try using the IntegralV3 factory instead.\n");
  fail();
  return 0;
  //  return new OneBodyDerivIntV3(this, bs1_, bs2_, &Int1eV3::hcore_1der);
}

Ref<TwoBodyInt>
IntegralLibint2::electron_repulsion()
{
  TwoBodyIntLibint2::ContractedGeminal bra, ket;  // these are dummy params for this evaluator anyway
  return new TwoBodyIntLibint2(this, bs1_, bs2_, bs3_, bs4_, storage_, erieval, bra, ket);
}

Ref<TwoBodyDerivInt>
IntegralLibint2::electron_repulsion_deriv()
{
  throw FeatureNotImplemented("IntegralLibint2::electron_repulsion_deriv() is not implemented yet", __FILE__, __LINE__);
}

Ref<TwoBodyInt>
IntegralLibint2::grt(const Ref<IntParams>& p)
{
  throw FeatureNotImplemented("IntegralLibint2::grt() is not implemented yet", __FILE__, __LINE__);
  //return new TwoBodyIntLibint2(this, bs1_, bs2_, bs3_, bs4_, storage_, grteval);
}

Ref<TwoBodyInt>
IntegralLibint2::g12(const Ref<IntParams>& params)
{
  Ref<IntParamsG12> params_cast;
  params_cast << params;
  if (params_cast.null())
    throw ProgrammingError("IntegralLibint2::g12() -- type of params does not match callback",__FILE__,__LINE__);
  return new TwoBodyIntLibint2(this, bs1_, bs2_, bs3_, bs4_, storage_,
                               g12eval, params_cast->ket(), params_cast->ket());
}

void
IntegralLibint2::set_basis(const Ref<GaussianBasisSet> &b1,
                           const Ref<GaussianBasisSet> &b2,
                           const Ref<GaussianBasisSet> &b3,
                           const Ref<GaussianBasisSet> &b4)
{
  free_transforms();
  Integral::set_basis(b1,b2,b3,b4);
  check_fullgencon();
  initialize_transforms();
}

void
IntegralLibint2::free_transforms()
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
IntegralLibint2::initialize_transforms()
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
    st_ = new SphericalTransformLibint2**[maxl_+1];
    ist_ = new ISphericalTransformLibint2**[maxl_+1];;
    int i,j;
    for (i=0; i<=maxl_; i++) {
      st_[i] = new SphericalTransformLibint2*[i/2+1];
      ist_[i] = new ISphericalTransformLibint2*[i/2+1];
      for (j=0; j<=i/2; j++) {
        st_[i][j] = new SphericalTransformLibint2(i,i-2*j);
        ist_[i][j] = new ISphericalTransformLibint2(i,i-2*j);
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
IntegralLibint2::check_fullgencon() const
{
  if ( has_fullgencon(bs1_) ||
       has_fullgencon(bs2_) ||
       has_fullgencon(bs3_) ||
       has_fullgencon(bs4_) ) {
      throw std::runtime_error("IntegralLibint2 cannot handle basis sets with fully general contractions yet, try IntegralV3 instead");
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
