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
    /// Reports symmetry with respect to the permutation of function in bra with function in ket for particle \c p
    /// +1 corresponds to hermitian, -1 to antihermitian, 0  to nonhermitian
    virtual int perm_symm(unsigned int p) const =0;
    /// Reports symmetry with respect to permutation of particles i and j in bra, or ket
    virtual int perm_symm(unsigned int i, unsigned int j) const =0;
  };

  /** Describes permutational properties (hermiticity) of one-body operators */
  class OneBodyOperDescr : public OperatorDescr {
  public:
    OneBodyOperDescr(int perm);

    /// Implementation of OperatorDescr::num_particles()
    unsigned int num_particles() const;
    /// Implementation of OperatorDescr::perm_symm()
    int perm_symm(unsigned int i) const;
    /// Implementation of OperatorDescr::perm_symm()
    int perm_symm(unsigned int i, unsigned int j) const;

  private:
    int perm_;
  };

  /** Describes permutational properties (hermiticity, Bose/Fermi) of a two-body operator */
  class TwoBodyOperDescr : public OperatorDescr {
  public:
    TwoBodyOperDescr(int perm_p1, int perm_p2, int perm_p12);

    /// Implementation of OperatorDescr::num_particles()
    unsigned int num_particles() const;
    /// Implementation of OperatorDescr::perm_symm()
    int perm_symm(unsigned int i) const;
    /// Implementation of OperatorDescr::perm_symm()
    int perm_symm(unsigned int i, unsigned int j) const;

  private:
    int perm_p1_;
    int perm_p2_;
    int perm_p12_;

  };

  /** Describes one-body operators.
   */
  struct OneBodyOper {
      /**
       * Types of one-body operators, includes various context-dependent "projectors", such as 1-RDM, etc.
       * "True" operators have type >= 0.
       */
    enum type {
      gamma = -1, //!< 1-body reduced density matrix
      T = 0,      //!< (nonrelativitic) kinetic energy
      V = 1,      //!< nuclear (Coulomb) potential
      h = 2,      //!< core Hamiltonian = T+V
      J = 3,      //!< (electronic) Coulomb
      K = 4,      //!< (electronic) exchange
      F = 5,      //!< Fock operator
      hJ = 6,     //!< h+J
      mu_x = 7,   //!< x component of electric dipole moment
      mu_y = 8,   //!< y component of electric dipole moment
      mu_z = 9,   //!< z component of electric dipole moment
      q_xx = 10,   //!< xx component of quadrupole moment
      q_xy = 11,  //!< xy component of quadrupole moment
      q_xz = 12,  //!< xz component of quadrupole moment
      q_yy = 13,  //!< yy component of quadrupole moment
      q_yz = 14,  //!< yz component of quadrupole moment
      q_zz = 15,  //!< zz component of quadrupole moment
      pVp = 16,   //!< \$f \underline{\hat{p}} \cdot V \underline{\hat{p}} \f$
      pxVp_x = 17,//!< x component of \$f \underline{\hat{p}} \cross V \underline{\hat{p}} \f$
      pxVp_y = 18,//!< y component of \$f \underline{\hat{p}} \cross V \underline{\hat{p}} \f$
      pxVp_z = 19,//!< z component of \$f \underline{\hat{p}} \cross V \underline{\hat{p}} \f$
      p4 = 20,      //!< \f$ (\underline{\hat{p}} \cdot \underline{\hat{p}})^2 \f$
      Nabla_x = 21,   //!< x component of Nabla operator ( \f$ \equiv i \hat{p}_x \f$ )
      Nabla_y = 22,   //!< y component of Nabla operator ( \f$ \equiv i \hat{p}_y \f$ )
      Nabla_z = 23,   //!< z component of Nabla operator ( \f$ \equiv i \hat{p}_z \f$ )
      iL_x = 24,   //!< x component of negative imaginary part of angular momentum ( \f$ \equiv i \hat{L}_x \f$ )
      iL_y = 25,   //!< y component of negative imaginary part of angular momentum ( \f$ \equiv i \hat{L}_y \f$ )
      iL_z = 26,    //!< z component of negative imaginary part of angular momentum ( \f$ \equiv i \hat{L}_z \f$ )
      invalid = 27
    };

    /// The max number of such types
    static int max_ntypes;

    /// Returns a descriptor for integral type t
    static Ref<OneBodyOperDescr> descr(OneBodyOper::type t);

    /// converts type to string
    static std::string to_string(type t);
  };

  /// Known one-body operator sets
  struct OneBodyOperSet {
      /**
       * one-body operator sets (\sa OneBodyOper::type)
       */
    enum type {
      T,  //!< {T}
      V,  //!< {V}
      h,  //!< {h}
      mu, //!< {mu_x, mu_y, mu_z}
      q,  //!< {q_xx, q_xy, q_xz, q_yy, q_yz, q_zz}
      pVp,//!< {pVp}
      p4  //!< {p4}
    };
  };

  /// Describes sets of two-body operators (\sa OneBodyOper)
  template <OneBodyOperSet::type Type> struct OneBodyOperSetTypeMap;
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::T> {
    static const int size = 1;
    static OneBodyOper::type value[];
  };
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::V> {
    static const int size = 1;
    static OneBodyOper::type value[];
  };
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::h> {
    static const int size = 1;
    static OneBodyOper::type value[];
  };
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::mu> {
    static const int size = 3;
    static OneBodyOper::type value[];
  };
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::q> {
    static const int size = 6;
    static OneBodyOper::type value[];
  };
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::pVp> {
    static const int size = 1;
    static OneBodyOper::type value[];
  };
  template <> struct OneBodyOperSetTypeMap<OneBodyOperSet::p4> {
    static const int size = 1;
    static OneBodyOper::type value[];
  };

  /// runtime version of OneBodyOperSetTypeMap
  class OneBodyOperSetDescr : public RefCount {
    public:
      static Ref<OneBodyOperSetDescr> instance(OneBodyOperSet::type oset);
      int size() const { return size_; }
      OneBodyOper::type opertype(unsigned int o) const;
      unsigned int opertype(OneBodyOper::type o) const;
    private:
      OneBodyOperSetDescr(int size,
                          const OneBodyOper::type* value) :
                          size_(size), value_(value) { }

      int size_;
      const OneBodyOper::type* value_;
  };

  /// which parameter set needed to specify the operator set?
  template <OneBodyOperSet::type Type> struct OneBodyIntParamsType;
  template <> struct OneBodyIntParamsType<OneBodyOperSet::T> {
    typedef IntParamsVoid value;
  };
  template <> struct OneBodyIntParamsType<OneBodyOperSet::V> {
    typedef IntParamsVoid value;
  };
  template <> struct OneBodyIntParamsType<OneBodyOperSet::h> {
    typedef IntParamsVoid value;
  };
  template <> struct OneBodyIntParamsType<OneBodyOperSet::mu> {
    typedef IntParamsOrigin value;
  };
  template <> struct OneBodyIntParamsType<OneBodyOperSet::q> {
    typedef IntParamsOrigin value;
  };
  template <> struct OneBodyIntParamsType<OneBodyOperSet::pVp> {
    typedef IntParamsVoid value;
  };
  template <> struct OneBodyIntParamsType<OneBodyOperSet::p4> {
    typedef IntParamsVoid value;
  };

  ////////////////////////////////////////////////////////////
  /// Two-Body Operators
  ////////////////////////////////////////////////////////////

  /** Describes types of two-body integrals
    */
  struct TwoBodyOper {
      /**
       * types of known two-body operators
       */
    enum type {
      eri =0,                  //!< two-body Coulomb \f$ 1/r_{12} \f$
      r12 =1,                  //!< interelectronic distance \f$ r_{12} \f$
      r12t1 =2,                //!< \f$ [r_{12},\hat{T}_1] \f$
      r12t2 =3,                //!< \f$ [r_{12},\hat{T}_2] \f$
      r12_0_g12 =4,            //!< (contracted) Gaussian geminal, \f$ g_{12} \equiv \sum_i c_i \exp[- \alpha_i r_{12}^2 ]\f$
      r12_m1_g12 =5,           //!< (contracted) Gaussian geminal over Coulomb,  \f$ r_{12}^{-1} g_{12} \f$
      t1g12 =6,                //!< \f$ [\hat{T}_1, g_{12}] \f$
      t2g12 =7,                //!< \f$ [\hat{T}_2, g_{12}] \f$
      g12t1g12 =8,             //!< \f$ [ g_{12}, [\hat{T}_1, g_{12}] ] \f$
      g12p4g12_m_g12t1g12t1 =9,//!< \f$ [g_{12}, [\hat{p}^4_1 + \hat{p}^4_2, g_{12}]] - 2 [g_{12}, [\hat{T}_1 + \hat{T}_2, g_{12}]](\hat{T}_1 + \hat{T}_2) \f$
      anti_g12g12 =10,         //!< anti_g12g12
      delta =11                //!< \f$ \delta_3({\bf r}_1 - {\bf r}_2) \f$
    };
    /// The max number of such types
    static int max_ntypes;
    /// Returns a descriptor for integral type t
    static Ref<TwoBodyOperDescr> descr(TwoBodyOper::type t);
    /// converts type to string representation
    static std::string to_string(type t);
    /// converts string representation to type
    static type to_type(const std::string& key);
  };

  /// Known two-body operator sets
  struct TwoBodyOperSet {
    enum type {
      ERI,       //!< {eri}
      R12,       //!< {eri, r12, r12t1, r12t2}
      G12,       //!< {eri, r12_0_g12, r12_m1_g12, t1g12, t2g12, g12t1g12}
      G12NC,     //!< {eri, r12_0_g12, r12_m1_g12, g12t1g12, anti_g12g12}
      G12DKH,    //!< {g12p4g12_m_g12t1g12t1}
      R12_0_G12, //!< {r12_0_g12}
      R12_m1_G12,//!< {r12_m1_g12}
      G12_T1_G12,//!< {g12t1g12}
      DeltaFunction //!< {delta}
    };

    /// converts type to string representation
    static std::string to_string(type t);
    /// converts string representation to type
    static type to_type(const std::string& key);

    /// maps TwoBodyOper::type to type
    /// @note only succeeds if the map is unequivocal
    static type to_type(TwoBodyOper::type oper);

  };

  /// Describes sets of two-body operators (\sa TwoBodyOper)
  template <TwoBodyOperSet::type Type> struct TwoBodyOperSetTypeMap;
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::ERI> {
    static const int size = 1;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::R12> {
    static const int size = 4;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::G12> {
    static const int size = 6;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::G12NC> {
    static const int size = 5;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::G12DKH> {
    static const int size = 1;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::R12_0_G12> {
    static const int size = 1;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::R12_m1_G12> {
    static const int size = 1;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::G12_T1_G12> {
    static const int size = 1;
    static TwoBodyOper::type value[];
    static std::string key;
  };
  template <> struct TwoBodyOperSetTypeMap<TwoBodyOperSet::DeltaFunction> {
    static const int size = 1;
    static TwoBodyOper::type value[];
    static std::string key;
  };

  /// runtime version of TwoBodyOperSetTypeMap
  class TwoBodyOperSetDescr : public RefCount {
    public:
      static Ref<TwoBodyOperSetDescr> instance(TwoBodyOperSet::type oset);
      int size() const { return size_; }
      TwoBodyOper::type opertype(unsigned int o) const;
      unsigned int opertype(TwoBodyOper::type o) const;
      std::string key() const { return key_; }
    private:
      TwoBodyOperSetDescr(int size,
                          const TwoBodyOper::type* value,
                          const std::string& key);
      int size_;
      const TwoBodyOper::type* value_;
      std::string key_;
  };

  /// which parameter set needed to specify the operator set?
  template <TwoBodyOperSet::type Type> struct TwoBodyIntParamsType;
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::ERI> {
    typedef IntParamsVoid value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::R12> {
    typedef IntParamsVoid value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::G12> {
    typedef IntParamsG12 value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::G12NC> {
    typedef IntParamsG12 value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::G12DKH> {
    typedef IntParamsG12 value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::R12_0_G12> {
    typedef IntParamsG12 value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::R12_m1_G12> {
    typedef IntParamsG12 value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::G12_T1_G12> {
    typedef IntParamsG12 value;
  };
  template <> struct TwoBodyIntParamsType<TwoBodyOperSet::DeltaFunction> {
    typedef IntParamsVoid value;
  };


}

#endif

