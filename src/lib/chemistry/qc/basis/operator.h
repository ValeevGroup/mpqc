//
// operator.h
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

#ifdef __GNUG__
#pragma interface
#endif

#include <util/ref/ref.h>
#include <math/distarray4/distarray4.h>
#include <chemistry/qc/basis/intparams.h>

#ifndef _chemistry_qc_basis_operator_h
#define _chemistry_qc_basis_operator_h

namespace sc {

  /** For an operator (e.g. <ij|1/r12|kl> )
      OperatorDescr gives basic properties with respect to permutations, etc.
  */
  class OperatorDescr : public RefCount {
  public:
    OperatorDescr() {}
    virtual ~OperatorDescr() {}

    virtual unsigned int num_particles() const =0;
    /// number of functions for particle i
    virtual unsigned int num_functions(unsigned int i) const =0;
    /// +1 corresponds to hermitian, -1 to antihermitian, 0  to nonhermitian
    virtual int perm_symm(unsigned int i) const =0;
    /// Symmetry with respect to permutation of particles i and j
    virtual int perm_symm(unsigned int i, unsigned int j) const =0;
  };

  /** Two-body operators with n functions per particle and given symmetry properties */
  class TwoBodyOperDescr : public OperatorDescr {
  public:
    TwoBodyOperDescr(unsigned int n, int perm_p1, int perm_p2, int perm_p12);

    /// Implementation of OperatorDescr::num_particles()
    unsigned int num_particles() const;
    /// Implementation of OperatorDescr::num_functions()
    unsigned int num_functions(unsigned int i) const;
    /// Implementation of OperatorDescr::perm_symm()
    int perm_symm(unsigned int i) const;
    /// Implementation of OperatorDescr::perm_symm()
    int perm_symm(unsigned int i, unsigned int j) const;

  private:
    unsigned int n_;
    int perm_p1_;
    int perm_p2_;
    int perm_p12_;

  };

  /** Types of two-body integrals that TwoBodyInt understands:
      eri -- integral of \f$r_{12}^{-1}\f$,
      r12 -- integral of \f$r_{12}\f$,
      r12ti -- integral of \f$[r_{12},\hat{T}_i]\f$,
      r12_0_g12 -- integral of \f$ g_{12}=\exp(-\gamma r_{12}^2) \f$,
      r12_m1_g12 -- integral of \f$g_{12}/r_{12}\f$,
      tig12 -- integral of \f$[\hat{T}_i,g_{12}]\f$,
      g12t1g12 -- integral of \f$[g_{12},[\hat{T}_1,g_{12}]]\f$,
      g12p4g12_m_g12t1g12t1 -- integral of
         \f$[g_{12}, [\hat{p}^4_1 + \hat{p}^4_2, g_{12}]] - 2 [g_{12}, [\hat{T}_1 + \hat{T}_2, g_{12}]](\hat{T}_1 + \hat{T}_2)\f$,
      anti_g12g12 -- integral of
       */
  struct TwoBodyOper {
    enum type { eri =0, r12 =1, r12t1 =2, r12t2 =3,
      r12_0_g12 =4, r12_m1_g12 =5, t1g12 =6, t2g12 =7,
      g12t1g12 =8, g12p4g12_m_g12t1g12t1 =9, anti_g12g12 =10,
      r12_0_gg12 =11, r12_m1_gg12 =12, gg12t1gg12 =13};
    /// The max number of such types
    static const int max_ntypes = DistArray4::max_num_te_types;
    /// Returns a descriptor for integral type t
    static const Ref<TwoBodyOperDescr>& descr(TwoBodyOper::type t);
  };

  /// Known operator sets
  struct TwoBodyOperSet {
    enum type {
      ERI,
      R12,
      G12,
      G12NC,
      GenG12,
      G12DKH
    };
  };

  /// Describes operator sets
  template <TwoBodyOperSet::type Type> struct OperSetTypeMap;
  template <> struct OperSetTypeMap<TwoBodyOperSet::ERI> {
    static const int size = 1;
    static TwoBodyOper::type value[];
  };
  template <> struct OperSetTypeMap<TwoBodyOperSet::R12> {
    static const int size = 4;
    static TwoBodyOper::type value[];
  };
  template <> struct OperSetTypeMap<TwoBodyOperSet::G12> {
    static const int size = 6;
    static TwoBodyOper::type value[];
  };
  template <> struct OperSetTypeMap<TwoBodyOperSet::G12NC> {
    static const int size = 5;
    static TwoBodyOper::type value[];
  };
  template <> struct OperSetTypeMap<TwoBodyOperSet::GenG12> {
    static const int size = 4;
    static TwoBodyOper::type value[];
  };
  template <> struct OperSetTypeMap<TwoBodyOperSet::G12DKH> {
    static const int size = 1;
    static TwoBodyOper::type value[];
  };

  /// runtime version of OperSetTypeMap
  class TwoBodyOperSetDescr : public RefCount {
    public:
      static Ref<TwoBodyOperSetDescr> instance(TwoBodyOperSet::type oset);
      int size() const { return size_; }
      TwoBodyOper::type opertype(unsigned int o) const;
      unsigned int opertype(TwoBodyOper::type o) const;
    private:
      TwoBodyOperSetDescr(int size,
                          const TwoBodyOper::type* value);
      int size_;
      const TwoBodyOper::type* value_;
  };

  /// which parameter set needed to specify the operator set?
  template <TwoBodyOperSet::type Type> struct IntParamsType;
  template <> struct IntParamsType<TwoBodyOperSet::ERI> {
    typedef IntParamsVoid value;
  };
  template <> struct IntParamsType<TwoBodyOperSet::R12> {
    typedef IntParamsVoid value;
  };
  template <> struct IntParamsType<TwoBodyOperSet::G12> {
    typedef IntParamsG12 value;
  };
  template <> struct IntParamsType<TwoBodyOperSet::G12NC> {
    typedef IntParamsG12 value;
  };
  template <> struct IntParamsType<TwoBodyOperSet::GenG12> {
    typedef IntParamsGenG12 value;
  };
  template <> struct IntParamsType<TwoBodyOperSet::G12DKH> {
    typedef IntParamsG12 value;
  };


}

#endif

