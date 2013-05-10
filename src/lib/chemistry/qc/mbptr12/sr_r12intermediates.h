//
// singlereference_r12_intermediates.h
//
// Copyright (C) 2013 Edward Valeev
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
#pragma implementation
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediates_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediates_h

#if defined(HAVE_MPQC3_RUNTIME)
# include <tiled_array.h>
#else
# error "sr_r12intermediates.h requires MPQC3 runtime, but it is not available"
#endif

#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <../bin/mpqc/mpqcinit.h>

namespace TA = TiledArray;

namespace sc {

  /// Tile of a 4-index tensor that's "evaluated" when needed by reading from DistArray4
  template <typename T>
  class DA4_Tile {
  private:
    TA::Array<T, 4, DA4_Tile>* owner_;
    std::array<std::size_t, 4> index_;
    Ref<DistArray4> darray4_;
    int te_type_;

  public:
    typedef T value_type;
    typedef TA::Tensor<T> eval_type;
    typedef typename eval_type::range_type range_type;

    DA4_Tile() { }

    DA4_Tile(TA::Array<T, 4, DA4_Tile>* owner,
             const std::array<std::size_t, 4>& index,
             const Ref<DistArray4>& darray4,
             int te_type) :
      owner_(owner), index_(index), darray4_(darray4), te_type_(te_type)
    { }

    operator TA::Tensor<T> () const;

    /// \tparam Archive The serialization archive type
    /// \param ar The serialization archive
    template <typename Archive>
    void serialize(Archive& ar) {
      assert(false);
    }

  };

  /// Tile of a <34> slice of <1234> that's "evaluated" when needed by reading from DistArray4 holding pqrs.
  template <typename T>
  class DA4_Tile34 {
  private:
    TA::Array<TA::Tensor<T>, 2, DA4_Tile34 >* owner_;
    std::array<std::size_t, 2> index_;
    Ref<DistArray4> darray4_;
    int te_type_;

  public:
    typedef T value_type;
    typedef TA::Tensor<TA::Tensor<T> > eval_type;
    typedef typename eval_type::range_type range_type;
    typedef T numeric_type;

    DA4_Tile34() { }

    DA4_Tile34(TA::Array<TA::Tensor<T>, 2, DA4_Tile34 >* owner,
               const std::array<std::size_t, 2>& index,
               const Ref<DistArray4>& darray4,
               int te_type) :
      owner_(owner), index_(index), darray4_(darray4), te_type_(te_type)
    { }

    /// conversion to eval_type
    operator TA::Tensor<TA::Tensor<T> > () const;

    /// \tparam Archive The serialization archive type
    /// \param ar The serialization archive
    template <typename Archive>
    void serialize(Archive& ar) {
      assert(false);
    }

  };

  namespace expressions {
    template <typename ArgType, bool Transpose> struct trace_tensor2_op;
    template <typename ArgType, bool Transpose> struct diag_tensor2_op;
  };


  /// SingleReference_R12Intermediates computes R12/F12 intermediates using MPQC3 runtime
  /// @tparam T the numeric type supporting all tensors
  template <typename T>
  class SingleReference_R12Intermediates {
    public:

      typedef TA::Array<T, 4 > TArray4; // Tile = Tensor<T>
      typedef TA::Array<T, 4, DA4_Tile<T> > TArray4d; // Tile = DA4_Tile<T>
      typedef TA::Array<T, 2> TArray2; // Tile = Tensor<T>
      typedef TA::Array<TA::Tensor<T>, 2> TArray22; // Tile = Tensor<Tensor<T>>
      typedef TA::Array<TA::Tensor<T>, 2, DA4_Tile34<T> > TArray22d; // Tile = Tensor<Tensor<T>>

      /**
       * Constructs an SingleReference_R12Intermediates object. There can be many.
       *
       * @param world MADNESS world in which the computed Arrays will live
       * @param r12world R12WavefunctionWorld that provides computation of integrals and OrbitalSpace objects
       */
      SingleReference_R12Intermediates(madness::World& world,
                                       const Ref<R12WavefunctionWorld>& r12world) :
                                         world_(world),
                                         r12world_(r12world)
      {
      }
      ~SingleReference_R12Intermediates() {
        r12world_ = 0;
      }

      /** computes diagonal (spin-restricted, for now) V intermediate
      * @return \f$ V^{ij}_{ij} \f$ and \f$ V^{ij}_{ji} \f$, respectively
      */
      std::pair<TArray2,TArray2> V_diag();

      /** computes diagonal (spin-restricted, for now) X intermediate
      * @return \f$ X^{ij}_{ij} \f$ and \f$ X^{ij}_{ji} \f$, respectively
      */
      std::pair<TArray2,TArray2> X_diag();

      /** computes diagonal (spin-restricted, for now) B intermediate
      * @return \f$ B^{ij}_{ij} \f$ and \f$ B^{ij}_{ji} \f$, respectively
      */
      std::pair<TArray2,TArray2> B_diag();

      /** computes 1-particle density intermediate
      * @return \f$ \gamma^{p}_{q} \f$, respectively
      */
      TArray2 rdm1();

      /**
       * provides T1 amplitude tensor
       * @param t1 act_occ by act_vir matrix
       */
      void set_T1(const RefSCMatrix& t1) { t1_ = t1; }
      /**
       * provides T1 CABS amplitude tensor
       * @param t1 act_occ by allvir (or CABS) matrix
       */
      void set_T1_cabs(const RefSCMatrix& t1_cabs) { t1_cabs_ = t1_cabs; }

      /// see _4()
      TArray4d& ijxy(const std::string& key);
      ///
      TArray22d& ij_xy(const std::string& key);
      /// see _2()
      TArray2& xy(const std::string& key);

      /** Creates a rank-4 Array, key describes the integral. A higher-level version of ijxy().
       *  key is similar to keys understood by ParsedTwoBodyInt and used by TwoBodyMOIntsRuntime,
       *  but with te_type embedded into key.
        * The following te_types are understood:
        *   - g : \f$ r_{12}^{-1} \f$
        *   - r : \f$ f(r_{12}) \f$
        *   - gr : \f$ r_{12}^{-1} f(r_{12}) \f$
        *   - rTr : \f$ [f(r_{12}), [ \hat{T}_1 , f(r_{12})]] \f$
        *
        * Here's the key index dictionary that can be used (\sa to_space):
        *   - i,j,k,l -> act_occ
        *   - m,n -> occ
        *   - a,b,c,d -> act_vir
        *   - e,f -> vir
        *   - p,q,r,s -> obs
        *   - p',q',r',s' -> cbs
        *   - a', b', c', d' -> cabs = cbs - obs
        *   - A', B', C', D' -> vir + cabs = cbs - occ
        *
        * Indices can also include a 1-particle operator transform (see Example2 and ParsedTransformedOrbitalSpaceKey)
        *
        * Example1: <i j|g|p a'> will create an array of <act_occ act_occ| r_{12}^{-1} |obs cabs> integrals
        * Example2: <i j_F(p')|g|p a'> will create an array of <act_occ cbs| r_{12}^{-1} |obs cabs> integrals
        * with the second index transformed using Fock matrix between act_occ and cbs spaces.
        *
        * Note that the indices are used to create the resulting tensor expression, hence these keys can be used
        * in composing expressions. For example, the full (non-diagonal) V intermediate without the CABS terms
        * can be computed as follows:
        * TArray4d V_ij_kl = _4("<i j|gr|k l>") - _4("<i j|g|p q>") * _4("<k l|r|p q>");
        *
        */
      TA::expressions::TensorExpression<TA::Tensor<T> > _4(const std::string& key);

      TA::expressions::TensorExpression<TA::Tensor<T> > _2(const std::string& key);

      //TA::expressions::TensorExpression<TA::Tensor< TA::Tensor<T> > > _22(const std::string& key);

    private:
      madness::World& world_;
      Ref<R12WavefunctionWorld> r12world_;

      // extra data
      RefSCMatrix t1_;
      RefSCMatrix t1_cabs_;

      // utilities

      /// computes a DistArray4 and converts into a TArray22
      std::shared_ptr<TArray22> ab_O_cd(const std::string& key,
                                        const int te_type);

      std::map<std::string, std::shared_ptr<TArray22> > tarray22_registry_;
      std::map<std::string, std::shared_ptr<TArray4d> > tarray4_registry_;
      std::map<std::string, std::shared_ptr<TArray2> > tarray2_registry_;
      std::map<std::string, std::shared_ptr<TArray22d> > tarray22d_registry_;

      // deprecated
      TArray22& _(const std::string& key);

      /**
       * converts index label to space label using the canonical dictionary documented in _()
       * @param index label
       * @return space label
       */
      std::string to_space(const std::string& index);

      typedef enum {
        ij=0, ji=1
      } ij_type;
      /// takes ij|o|ij or ij|o|ji
      template <typename Array22>
      TA::expressions::TensorExpression<TA::Tensor<T> > take(const Array22& ij_o_pq,
                                                             ij_type IJ) {
        //typedef typename Array22::value_type value_type;
          typedef TA::Tensor<TA::Tensor<T> > value_type;
            if (IJ == ij) {
              sc::expressions::diag_tensor2_op<value_type,false> diag_op;
              return make_unary_tensor(ij_o_pq("i,j"),
                                       diag_op);
            }
            // else: IJ == ji
            sc::expressions::diag_tensor2_op<value_type,true> diag_op;
            return make_unary_tensor(ij_o_pq("i,j"),
                                     diag_op);
          }

      /// computes ij|o1|pq . ij|o2|pq
      template <typename Array22>
      TA::expressions::TensorExpression<TA::Tensor<T> > dotket(const Array22& ij_o1_pq,
                                                               const Array22& ij_o2_pq,
                                                               bool transpose_o2_ket = false) {
        //typedef typename Array22::value_type value_type;
        typedef TA::Tensor<TA::Tensor<T> > value_type;
        if (transpose_o2_ket == false) {
          sc::expressions::trace_tensor2_op<value_type,false> trace_op;
          return make_binary_tensor(ij_o1_pq("i,j"), ij_o2_pq("i,j"),
                                    trace_op);
        }
        // else: transpose_o2_ket = true
        sc::expressions::trace_tensor2_op<value_type,false> trace_op;
        return make_binary_tensor(ij_o1_pq("i,j"), ij_o2_pq("j,i"),
                                  trace_op);
      }

  };

}; // end of namespace sc

#include <chemistry/qc/mbptr12/sr_r12intermediates_util.h>
#include <chemistry/qc/mbptr12/sr_r12intermediates_VXB_diag.h>

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
