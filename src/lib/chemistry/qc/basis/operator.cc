//
// operator.cc
//
// Copyright (C) 2007 Edward Valeev
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

#include <cassert>
#include <util/class/scexception.h>
#include <chemistry/qc/basis/operator.h>

using namespace sc;

////

int TwoBodyOper::max_ntypes = DistArray4::max_num_te_types;

TwoBodyOperDescr::TwoBodyOperDescr(unsigned int n, int perm_p1, int perm_p2, int perm_p12) :
  n_(n), perm_p1_(perm_p1), perm_p2_(perm_p2), perm_p12_(perm_p12)
{
}

unsigned int TwoBodyOperDescr::num_particles() const { return 2; }
unsigned int TwoBodyOperDescr::num_functions(unsigned int i) const { return n_; }
int TwoBodyOperDescr::perm_symm(unsigned int i) const {
  if (i == 1) return perm_p1_;
  if (i == 2) return perm_p2_;
  throw ProgrammingError("IntBodyIntTypeDescr::perm_symm(i) -- i must be 1 or 2",__FILE__,__LINE__);
}
int TwoBodyOperDescr::perm_symm(unsigned int i, unsigned int j) const {
  if ( (i == 1 && j == 2) ||
       (i == 2 && j == 1) ) return perm_p12_;
  throw ProgrammingError("IntBodyIntTypeDescr::perm_symm(i,j) -- i,j must be 1,2",__FILE__,__LINE__);
}

Ref<TwoBodyOperDescr>
TwoBodyOper::descr(TwoBodyOper::type type)
{
    static Ref<TwoBodyOperDescr> symm_type = new TwoBodyOperDescr(2,+1,+1,+1);
    static Ref<TwoBodyOperDescr> t1r12_inttype = new TwoBodyOperDescr(2,-1,+1,0);
    static Ref<TwoBodyOperDescr> t2r12_inttype = new TwoBodyOperDescr(2,+1,-1,0);
    switch(type) {
    case TwoBodyOper::r12t1:
    case TwoBodyOper::t1g12:
    return t1r12_inttype;

    case TwoBodyOper::r12t2:
    case TwoBodyOper::t2g12:
    return t2r12_inttype;

    case TwoBodyOper::eri:
    case TwoBodyOper::r12:
    case TwoBodyOper::r12_0_g12:
    case TwoBodyOper::r12_m1_g12:
    case TwoBodyOper::g12t1g12:
    case TwoBodyOper::anti_g12g12:
    case TwoBodyOper::r12_0_gg12:
    case TwoBodyOper::r12_m1_gg12:
    case TwoBodyOper::gg12t1gg12:
    case TwoBodyOper::g12p4g12_m_g12t1g12t1:
    return symm_type;
    }
    throw ProgrammingError("TwoBodyOper::descr() -- incorrect type");
}

TwoBodyOper::type OperSetTypeMap<TwoBodyOperSet::ERI>::value[] = {TwoBodyOper::eri};
TwoBodyOper::type OperSetTypeMap<TwoBodyOperSet::R12>::value[] = {TwoBodyOper::eri,
                                                            TwoBodyOper::r12,
                                                            TwoBodyOper::r12t1,
                                                            TwoBodyOper::r12t2};
TwoBodyOper::type OperSetTypeMap<TwoBodyOperSet::G12>::value[] = {TwoBodyOper::eri,
                                                            TwoBodyOper::r12_0_g12,
                                                            TwoBodyOper::r12_m1_g12,
                                                            TwoBodyOper::t1g12,
                                                            TwoBodyOper::t2g12,
                                                            TwoBodyOper::g12t1g12};
TwoBodyOper::type OperSetTypeMap<TwoBodyOperSet::G12NC>::value[] = {TwoBodyOper::eri,
                                                              TwoBodyOper::r12_0_g12,
                                                              TwoBodyOper::r12_m1_g12,
                                                              TwoBodyOper::g12t1g12,
                                                              TwoBodyOper::anti_g12g12};
TwoBodyOper::type OperSetTypeMap<TwoBodyOperSet::G12DKH>::value[] = {TwoBodyOper::g12p4g12_m_g12t1g12t1};

TwoBodyOperSetDescr::TwoBodyOperSetDescr(int size,
                                         const TwoBodyOper::type* value) :
  size_(size), value_(value)
  {}

Ref<TwoBodyOperSetDescr>
TwoBodyOperSetDescr::instance(TwoBodyOperSet::type oset)
{
  switch (oset) {
    case TwoBodyOperSet::ERI:
      return new TwoBodyOperSetDescr(OperSetTypeMap<TwoBodyOperSet::ERI>::size,
                                     OperSetTypeMap<TwoBodyOperSet::ERI>::value);
      break;
    case TwoBodyOperSet::G12:
      return new TwoBodyOperSetDescr(OperSetTypeMap<TwoBodyOperSet::G12>::size,
                                     OperSetTypeMap<TwoBodyOperSet::G12>::value);
      break;
    case TwoBodyOperSet::G12NC:
      return new TwoBodyOperSetDescr(OperSetTypeMap<TwoBodyOperSet::G12NC>::size,
                                     OperSetTypeMap<TwoBodyOperSet::G12NC>::value);
      break;
    case TwoBodyOperSet::G12DKH:
      return new TwoBodyOperSetDescr(OperSetTypeMap<TwoBodyOperSet::G12DKH>::size,
                                     OperSetTypeMap<TwoBodyOperSet::G12DKH>::value);
      break;

    default:
      assert(false);
  }
  return Ref<TwoBodyOperSetDescr>(); // unreachable
}

TwoBodyOper::type
TwoBodyOperSetDescr::opertype(unsigned int o) const
{
  assert(o < size_);
  return value_[o];
}

unsigned int
TwoBodyOperSetDescr::opertype(TwoBodyOper::type o) const
{
  for(int i=0; i<size_; ++i)
    if (o == value_[i])
      return i;
  abort();  // should be unreachable
}
