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
#include <util/misc/scexception.h>
#include <chemistry/qc/basis/operator.h>

using namespace sc;

////

int OneBodyOper::max_ntypes = static_cast<int>(OneBodyOper::invalid);

Ref<OneBodyOperDescr>
OneBodyOper::descr(OneBodyOper::type t) {
  static Ref<OneBodyOperDescr> hermitian = new OneBodyOperDescr(+1);
  static Ref<OneBodyOperDescr> antihermitian = new OneBodyOperDescr(-1);
  switch(t) {
    case Nabla_x:
    case Nabla_y:
    case Nabla_z:
    case iL_x:
    case iL_y:
    case iL_z:
      return antihermitian;

    case gamma:
    case T:
    case V:
    case h:
    case J:
    case K:
    case F:
    case mu_x:
    case mu_y:
    case mu_z:
    case q_xx:
    case q_xy:
    case q_xz:
    case q_yy:
    case q_yz:
    case q_zz:
    case p4:
    case pVp:
    case pxVp_x:
    case pxVp_y:
    case pxVp_z:
    case phi:
    case dphi_x:
    case dphi_y:
    case dphi_z:
    case ddphi_xx:
    case ddphi_xy:
    case ddphi_xz:
    case ddphi_yy:
    case ddphi_yz:
    case ddphi_zz:
      return hermitian;

    default:
      throw ProgrammingError("TwoBodyOper::descr() -- incorrect type");
  }

  return Ref<OneBodyOperDescr>(); // shut up stupid compiler
}

std::string
OneBodyOper::to_string(OneBodyOper::type t) {
  switch (t) {
    case gamma: return "gamma"; break;
    case S: return "S"; break;
    case T: return "T"; break;
    case V: return "V"; break;
    case h: return "h"; break;
    case J: return "J"; break;
    case K: return "K"; break;
    case F: return "F"; break;
    case hJ: return "hJ"; break;
    case mu_x: return "mu_x"; break;
    case mu_y: return "mu_y"; break;
    case mu_z: return "mu_z"; break;
    case q_xx: return "q_xx"; break;
    case q_xy: return "q_xy"; break;
    case q_xz: return "q_xz"; break;
    case q_yy: return "q_yy"; break;
    case q_yz: return "q_yz"; break;
    case q_zz: return "q_zz"; break;
    case pVp: return "pVp"; break;
    case pxVp_x: return "pxVp_x"; break;
    case pxVp_y: return "pxVp_y"; break;
    case pxVp_z: return "pxVp_z"; break;
    case p4: return "p4"; break;
    case Nabla_x: return "Nabla_x"; break;
    case Nabla_y: return "Nabla_y"; break;
    case Nabla_z: return "Nabla_z"; break;
    case iL_x: return "iL_x"; break;
    case iL_y: return "iL_y"; break;
    case iL_z: return "iL_z"; break;
    case phi: return "phi"; break;
    case dphi_x: return "dphi_x"; break;
    case dphi_y: return "dphi_y"; break;
    case dphi_z: return "dphi_z"; break;
    case ddphi_xx: return "ddphi_xx"; break;
    case ddphi_xy: return "ddphi_xy"; break;
    case ddphi_xz: return "ddphi_xz"; break;
    case ddphi_yy: return "ddphi_yy"; break;
    case ddphi_yz: return "ddphi_yz"; break;
    case ddphi_zz: return "ddphi_zz"; break;
    default:
      std::ostringstream oss;
      oss << "OneBodyOper::to_string: unknown type " << t;
      throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
}

OneBodyOper::type
OneBodyOper::to_type(const std::string& key) {
  for(int t=0; t<max_ntypes; ++t) {
    type tt = static_cast<type>(t);
    if (key == to_string(tt))
      return tt;
  }
  MPQC_ASSERT(false); // should be unreachable
  return invalid;
}

std::string
OneBodyOperSet::to_string(type t) {
  switch (t) {
    case S: return "S"; break;
    case T: return "T"; break;
    case V: return "V"; break;
    case h: return "h"; break;
    case mu: return "mu"; break;
    case q: return "q"; break;
    case phi: return "phi"; break;
    case dphi: return "dphi"; break;
    case ddphi: return "ddphi"; break;
    case p4: return "p4"; break;
    case pVp: return "pVp"; break;
    default: MPQC_ASSERT(false);
  }
  MPQC_ASSERT(false);
  return "";
}

OneBodyOperSet::type
OneBodyOperSet::to_type(const std::string& key) {
  for(int t=0; t<invalid; ++t) {
      type tt = static_cast<type>(t);
      if (key == to_string(tt))
        return tt;
    }
    MPQC_ASSERT(false); // should be unreachable
    return invalid;
}

OneBodyOperSet::type
OneBodyOperSet::to_type(OneBodyOper::type oper) {
  switch (oper) {
    case OneBodyOper::S: return S;
    case OneBodyOper::T: return T;
    case OneBodyOper::V: return V;
    case OneBodyOper::h: return h;
    case OneBodyOper::mu_x:
    case OneBodyOper::mu_y:
    case OneBodyOper::mu_z:
      return mu;
    case OneBodyOper::q_xx:
    case OneBodyOper::q_xy:
    case OneBodyOper::q_xz:
    case OneBodyOper::q_yy:
    case OneBodyOper::q_yz:
    case OneBodyOper::q_zz:
      return q;
    case OneBodyOper::phi: return phi;
    case OneBodyOper::dphi_x:
    case OneBodyOper::dphi_y:
    case OneBodyOper::dphi_z:
      return dphi;
    case OneBodyOper::ddphi_xx:
    case OneBodyOper::ddphi_xy:
    case OneBodyOper::ddphi_xz:
    case OneBodyOper::ddphi_yy:
    case OneBodyOper::ddphi_yz:
    case OneBodyOper::ddphi_zz:
      return ddphi;
    case OneBodyOper::p4: return p4;
    case OneBodyOper::pVp: return pVp;
    default: MPQC_ASSERT(false);
  }
  return invalid;
}

OneBodyOperDescr::OneBodyOperDescr(int perm) : perm_(perm) {
}

unsigned int
OneBodyOperDescr::num_particles() const {
  return 1;
}

int
OneBodyOperDescr::perm_symm(unsigned int i) const {
  MPQC_ASSERT(i == 0);
  return perm_;
}

int
OneBodyOperDescr::perm_symm(unsigned int i, unsigned int j) const {
  throw ProgrammingError("OneBodyOperDescr::perm_symm(i,j) -- lacks meaning", __FILE__, __LINE__);
}

OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::S>::value[] = {OneBodyOper::S};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::T>::value[] = {OneBodyOper::T};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::V>::value[] = {OneBodyOper::V};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::h>::value[] = {OneBodyOper::h};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::mu>::value[] = {OneBodyOper::mu_x,
                                                                 OneBodyOper::mu_y,
                                                                 OneBodyOper::mu_z};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::q>::value[] = {OneBodyOper::q_xx,
                                                                OneBodyOper::q_xy,
                                                                OneBodyOper::q_xz,
                                                                OneBodyOper::q_yy,
                                                                OneBodyOper::q_yz,
                                                                OneBodyOper::q_zz};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::phi>::value[] = {OneBodyOper::phi};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::dphi>::value[] = {OneBodyOper::dphi_x,
                                                                          OneBodyOper::dphi_y,
                                                                          OneBodyOper::dphi_z};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::ddphi>::value[] = {OneBodyOper::ddphi_xx,
                                                                          OneBodyOper::ddphi_xy,
                                                                          OneBodyOper::ddphi_xz,
                                                                          OneBodyOper::ddphi_yy,
                                                                          OneBodyOper::ddphi_yz,
                                                                          OneBodyOper::ddphi_zz};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::pVp>::value[] = {OneBodyOper::pVp};
OneBodyOper::type OneBodyOperSetProperties<OneBodyOperSet::p4>::value[] = {OneBodyOper::p4};

std::string OneBodyOperSetProperties<OneBodyOperSet::S>::key("S");
std::string OneBodyOperSetProperties<OneBodyOperSet::T>::key("T");
std::string OneBodyOperSetProperties<OneBodyOperSet::V>::key("V");
std::string OneBodyOperSetProperties<OneBodyOperSet::h>::key("h");
std::string OneBodyOperSetProperties<OneBodyOperSet::mu>::key("mu");
std::string OneBodyOperSetProperties<OneBodyOperSet::q>::key("q");
std::string OneBodyOperSetProperties<OneBodyOperSet::phi>::key("phi");
std::string OneBodyOperSetProperties<OneBodyOperSet::dphi>::key("dphi");
std::string OneBodyOperSetProperties<OneBodyOperSet::ddphi>::key("ddphi");
std::string OneBodyOperSetProperties<OneBodyOperSet::pVp>::key("pVp");
std::string OneBodyOperSetProperties<OneBodyOperSet::p4>::key("p4");

Ref<OneBodyOperSetDescr>
OneBodyOperSetDescr::instance(OneBodyOperSet::type oset)
{
  switch (oset) {
    case OneBodyOperSet::S:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::S>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::S>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::S>::key);
      break;
    case OneBodyOperSet::T:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::T>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::T>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::S>::key);
      break;
    case OneBodyOperSet::V:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::V>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::V>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::V>::key);
      break;
    case OneBodyOperSet::h:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::h>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::h>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::h>::key);
      break;
    case OneBodyOperSet::mu:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::mu>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::mu>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::mu>::key);
      break;
    case OneBodyOperSet::q:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::q>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::q>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::q>::key);
      break;
    case OneBodyOperSet::phi:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::phi>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::phi>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::phi>::key);
      break;
    case OneBodyOperSet::dphi:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::dphi>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::dphi>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::dphi>::key);
      break;
    case OneBodyOperSet::ddphi:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::ddphi>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::ddphi>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::ddphi>::key);
      break;
    case OneBodyOperSet::pVp:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::pVp>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::pVp>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::pVp>::key);
      break;
    case OneBodyOperSet::p4:
      return new OneBodyOperSetDescr(OneBodyOperSetProperties<OneBodyOperSet::p4>::size,
                                     OneBodyOperSetProperties<OneBodyOperSet::p4>::value,
                                     OneBodyOperSetProperties<OneBodyOperSet::p4>::key);
      break;
    default:
      MPQC_ASSERT(false);
  }
  return Ref<OneBodyOperSetDescr>(); // unreachable
}

OneBodyOper::type
OneBodyOperSetDescr::opertype(unsigned int o) const
{
  MPQC_ASSERT(o < size_);
  return value_[o];
}

unsigned int
OneBodyOperSetDescr::opertype(OneBodyOper::type o) const
{
  for(int i=0; i<size_; ++i)
    if (o == value_[i])
      return i;
  abort();  // should be unreachable
}

////

int TwoBodyOper::max_ntypes = DistArray4::max_num_te_types;

TwoBodyOperDescr::TwoBodyOperDescr(int perm_p1, int perm_p2, int perm_p12) :
  perm_p1_(perm_p1), perm_p2_(perm_p2), perm_p12_(perm_p12)
{
}

unsigned int TwoBodyOperDescr::num_particles() const { return 2; }
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
    static Ref<TwoBodyOperDescr> symm_type = new TwoBodyOperDescr(+1,+1,+1);
    static Ref<TwoBodyOperDescr> t1r12_inttype = new TwoBodyOperDescr(-1,+1,0);
    static Ref<TwoBodyOperDescr> t2r12_inttype = new TwoBodyOperDescr(+1,-1,0);
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
    case TwoBodyOper::g12p4g12_m_g12t1g12t1:
    case TwoBodyOper::delta:
    return symm_type;
    }
    throw ProgrammingError("TwoBodyOper::descr() -- incorrect type");
}

std::string
TwoBodyOper::to_string(TwoBodyOper::type t) {
  switch (t) {
    case eri: return "g"; break;
    case r12: return "r12"; break;
    case r12t1: return "r12t1"; break;
    case r12t2: return "r12t2"; break;
    case r12_0_g12: return "g12"; break;
    case r12_m1_g12: return "g12g"; break;
    case g12t1g12: return "g12t1g12"; break;
    case t1g12: return "t1g12"; break;
    case t2g12: return "t2g12"; break;
    case delta: return "delta"; break;
    default:
      std::ostringstream oss;
      oss << "TwoBodyOper::to_string: unknown type " << t;
      throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
}

sc::TwoBodyOper::type sc::TwoBodyOper::to_type(const std::string& key) {
  if (key == "g")
    return TwoBodyOper::eri;
  if (key == "r12")
    return TwoBodyOper::r12;
  if (key == "r12t1")
    return TwoBodyOper::r12t1;
  if (key == "r12t2")
    return TwoBodyOper::r12t2;
  if (key == "g12")
    return TwoBodyOper::r12_0_g12;
  if (key == "g12g")
    return TwoBodyOper::r12_m1_g12;
  if (key == "g12t1g12")
    return TwoBodyOper::g12t1g12;
  if (key == "t1g12")
    return TwoBodyOper::t1g12;
  if (key == "t2g12")
    return TwoBodyOper::t2g12;
  if (key == "delta")
    return TwoBodyOper::delta;
  std::ostringstream oss;
  oss << "TwoBodyOper::to_type: unknown string " << key;
  throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
}

//////////////////////////////////////////////////////

TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::ERI>::value[] = {TwoBodyOper::eri};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::R12>::value[] = {TwoBodyOper::eri,
                                                            TwoBodyOper::r12,
                                                            TwoBodyOper::r12t1,
                                                            TwoBodyOper::r12t2};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::G12>::value[] = {TwoBodyOper::eri,
                                                            TwoBodyOper::r12_0_g12,
                                                            TwoBodyOper::r12_m1_g12,
                                                            TwoBodyOper::t1g12,
                                                            TwoBodyOper::t2g12,
                                                            TwoBodyOper::g12t1g12};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::G12NC>::value[] = {TwoBodyOper::eri,
                                                              TwoBodyOper::r12_0_g12,
                                                              TwoBodyOper::r12_m1_g12,
                                                              TwoBodyOper::g12t1g12,
                                                              TwoBodyOper::anti_g12g12};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::G12DKH>::value[] = {TwoBodyOper::g12p4g12_m_g12t1g12t1};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::R12_0_G12>::value[] = {TwoBodyOper::r12_0_g12};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::R12_m1_G12>::value[] = {TwoBodyOper::r12_m1_g12};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::G12_T1_G12>::value[] = {TwoBodyOper::g12t1g12};
TwoBodyOper::type TwoBodyOperSetProperties<TwoBodyOperSet::DeltaFunction>::value[] = {TwoBodyOper::delta};

std::string TwoBodyOperSetProperties<TwoBodyOperSet::ERI>::key("ERI");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::G12>::key("G12");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::G12NC>::key("G12'");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::G12DKH>::key("G12DKH");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::R12_0_G12>::key("R12_0_G12");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::R12_m1_G12>::key("R12_m1_G12");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::G12_T1_G12>::key("G12_T1_G12");
std::string TwoBodyOperSetProperties<TwoBodyOperSet::DeltaFunction>::key("Delta");

TwoBodyOperSetDescr::TwoBodyOperSetDescr(int size,
                                         const TwoBodyOper::type* value,
                                         const std::string& key) :
  size_(size), value_(value), key_(key)
  {}

Ref<TwoBodyOperSetDescr>
TwoBodyOperSetDescr::instance(TwoBodyOperSet::type oset)
{
  switch (oset) {
    case TwoBodyOperSet::ERI:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::ERI>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::ERI>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::ERI>::key);
      break;
    case TwoBodyOperSet::R12:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::G12>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12>::key);
      break;
    case TwoBodyOperSet::G12:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::G12>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12>::key);
      break;
    case TwoBodyOperSet::G12NC:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::G12NC>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12NC>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12NC>::key);
      break;
    case TwoBodyOperSet::G12DKH:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::G12DKH>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12DKH>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12DKH>::key);
      break;
    case TwoBodyOperSet::R12_0_G12:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::R12_0_G12>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::R12_0_G12>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::R12_0_G12>::key);
      break;
    case TwoBodyOperSet::R12_m1_G12:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::R12_m1_G12>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::R12_m1_G12>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::R12_m1_G12>::key);
      break;
    case TwoBodyOperSet::G12_T1_G12:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::G12_T1_G12>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12_T1_G12>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::G12_T1_G12>::key);
      break;
    case TwoBodyOperSet::DeltaFunction:
      return new TwoBodyOperSetDescr(TwoBodyOperSetProperties<TwoBodyOperSet::DeltaFunction>::size,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::DeltaFunction>::value,
                                     TwoBodyOperSetProperties<TwoBodyOperSet::DeltaFunction>::key);
      break;

    default:
      MPQC_ASSERT(false);
  }
  return Ref<TwoBodyOperSetDescr>(); // unreachable
}

TwoBodyOper::type
TwoBodyOperSetDescr::opertype(unsigned int o) const
{
  MPQC_ASSERT(o < size_);
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

std::string
sc::TwoBodyOperSet::to_string(type t) {
  Ref<TwoBodyOperSetDescr> descr = TwoBodyOperSetDescr::instance(t);
  return descr->key();
}

sc::TwoBodyOperSet::type
sc::TwoBodyOperSet::to_type(const std::string& key) {
  for(int s=int(TwoBodyOperSet::ERI);
      s<=int(TwoBodyOperSet::DeltaFunction);
      ++s) {
    if (key == to_string(static_cast<TwoBodyOperSet::type>(s)))
      return static_cast<TwoBodyOperSet::type>(s);
  }
  std::ostringstream oss;
  oss << "TwoBodyOperSet::to_type: unknown string " << key;
  throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
}

sc::TwoBodyOperSet::type
sc::TwoBodyOperSet::to_type(TwoBodyOper::type oper) {
  if (oper == TwoBodyOper::eri)
    return ERI;
  if (oper == TwoBodyOper::r12_0_g12)
    return R12_0_G12;
  if (oper == TwoBodyOper::r12_m1_g12)
    return R12_m1_G12;
  if (oper == TwoBodyOper::g12t1g12)
    return G12_T1_G12;
  if (oper == TwoBodyOper::delta)
    return DeltaFunction;

  std::ostringstream oss;
  oss << "TwoBodyOperSet::to_type: map from " << TwoBodyOper::to_string(oper) << " not defined";
  throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
}
