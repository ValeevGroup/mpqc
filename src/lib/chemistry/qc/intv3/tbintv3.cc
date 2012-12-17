//
// tbintv3.cc
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

#include <stdexcept>

#include <util/misc/scexception.h>
#include <chemistry/qc/intv3/tbintv3.h>
#include <chemistry/qc/basis/integral.h>

using namespace sc;

TwoBodyIntV3::TwoBodyIntV3(Integral*integral,
                           const Ref<GaussianBasisSet>& b1,
                           const Ref<GaussianBasisSet>& b2,
                           const Ref<GaussianBasisSet>& b3,
                           const Ref<GaussianBasisSet>& b4,
                           size_t storage):
  TwoBodyInt(integral,b1,b2,b3,b4),
  descr_(TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI))
{
  int2ev3_ = new Int2eV3(integral,b1,b2,b3,b4,0,storage);
  buffer_ = int2ev3_->buffer();
  integral_->adjust_storage(int2ev3_->used_storage());
}

TwoBodyIntV3::~TwoBodyIntV3()
{
  integral_->adjust_storage(-int2ev3_->used_storage());
}

void
TwoBodyIntV3::compute_shell(int is, int js, int ks, int ls)
{
  int2ev3_->set_redundant(redundant());
  int2ev3_->erep(is,js,ks,ls);
}

int
TwoBodyIntV3::log2_shell_bound(int is, int js, int ks, int ls)
{
  return int2ev3_->erep_4bound(is,js,ks,ls);
}

void
TwoBodyIntV3::set_integral_storage(size_t storage)
{
  int2ev3_->init_storage(storage);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyThreeCenterIntV3::TwoBodyThreeCenterIntV3(
    Integral*integral,
    const Ref<GaussianBasisSet>& b1,
    const Ref<GaussianBasisSet>& b2,
    const Ref<GaussianBasisSet>& b3,
    size_t storage):
  TwoBodyThreeCenterInt(integral,b1,b2,b3),
  descr_(TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI))
{
  Ref<GaussianBasisSet> null;
  int2ev3_ = new Int2eV3(integral,b1,b2,b3,null,0,storage);
  buffer_ = int2ev3_->buffer();
  integral_->adjust_storage(int2ev3_->used_storage());
}

TwoBodyThreeCenterIntV3::~TwoBodyThreeCenterIntV3()
{
  integral_->adjust_storage(-int2ev3_->used_storage());
}

void
TwoBodyThreeCenterIntV3::compute_shell(int is, int js, int ks)
{
  int2ev3_->set_redundant(redundant());
  int2ev3_->erep_3center(is,js,ks);
}

int
TwoBodyThreeCenterIntV3::log2_shell_bound(int is, int js, int ks)
{
  throw std::runtime_error("TwoBodyThreeCenterIntv3: doesn't support bounds");
  return 0;
}

void
TwoBodyThreeCenterIntV3::set_integral_storage(size_t storage)
{
  int2ev3_->init_storage(storage);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyTwoCenterIntV3::TwoBodyTwoCenterIntV3(
    Integral*integral,
    const Ref<GaussianBasisSet>& b1,
    const Ref<GaussianBasisSet>& b2,
    size_t storage):
  TwoBodyTwoCenterInt(integral,b1,b2),
  descr_(TwoBodyOperSetDescr::instance(TwoBodyOperSet::ERI))
{
  Ref<GaussianBasisSet> null;
  int2ev3_ = new Int2eV3(integral,b1,null,b2,null,0,storage);
  buffer_ = int2ev3_->buffer();
  integral_->adjust_storage(int2ev3_->used_storage());
}

TwoBodyTwoCenterIntV3::~TwoBodyTwoCenterIntV3()
{
  integral_->adjust_storage(-int2ev3_->used_storage());
}

void
TwoBodyTwoCenterIntV3::compute_shell(int is, int js)
{
  int2ev3_->set_redundant(redundant());
  int2ev3_->erep_2center(is,js);
}

int
TwoBodyTwoCenterIntV3::log2_shell_bound(int is, int js)
{
  throw std::runtime_error("TwoBodyTwoCenterIntv3: doesn't support bounds");
  return 0;
}

void
TwoBodyTwoCenterIntV3::set_integral_storage(size_t storage)
{
  int2ev3_->init_storage(storage);
}

bool
TwoBodyTwoCenterIntV3::cloneable()
{
  return true;
}

Ref<TwoBodyTwoCenterInt>
TwoBodyTwoCenterIntV3::clone()
{
  const size_t storage = int2ev3_->used_storage();
  return new TwoBodyTwoCenterIntV3(integral_, bs1_, bs2_, storage);
}

//////////////////////////////////////////////////////////////////////////

TwoBodyDerivIntV3::TwoBodyDerivIntV3(Integral*integral,
                                     const Ref<GaussianBasisSet>& b1,
                                     const Ref<GaussianBasisSet>& b2,
                                     const Ref<GaussianBasisSet>& b3,
                                     const Ref<GaussianBasisSet>& b4,
                                     size_t storage):
  TwoBodyDerivInt(integral,b1,b2,b3,b4)
{
  int2ev3_ = new Int2eV3(integral,b1,b2,b3,b4,1,storage);
  buffer_ = int2ev3_->buffer();
  integral_->adjust_storage(int2ev3_->used_storage());
}

TwoBodyDerivIntV3::~TwoBodyDerivIntV3()
{
  integral_->adjust_storage(-int2ev3_->used_storage());
}

void
TwoBodyDerivIntV3::compute_shell(int is, int js, int ks, int ls,
                                 DerivCenters&c)
{
  int center;
  der_centersv3_t dercenters;
  int sh[4], sz[4];

  sh[0]=is; sh[1]=js; sh[2]=ks; sh[3]=ls;

  int2ev3_->erep_all1der(sh,sz,&dercenters);

  c.clear();
  for (int i=0; i<dercenters.n; i++) {
      if (dercenters.cs[i] == int2ev3_->pcs1()) center = 0;
      else if (dercenters.cs[i] == int2ev3_->pcs2()) center = 1;
      else if (dercenters.cs[i] == int2ev3_->pcs3()) center = 2;
      else center = 3;
      c.add_center(center,dercenters.num[i]);
    }
  if (dercenters.n) {
      if (dercenters.ocs == int2ev3_->pcs1()) center = 0;
      else if (dercenters.ocs == int2ev3_->pcs2()) center = 1;
      else if (dercenters.ocs == int2ev3_->pcs3()) center = 2;
      else center = 3;
      c.add_omitted(center,dercenters.onum);
    }
}


int
TwoBodyDerivIntV3::log2_shell_bound(int is, int js, int ks, int ls)
{
  return int2ev3_->erep_4bound_1der(is,js,ks,ls);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
