//
// singles_casscf.h
//
// Copyright (C) 2014 Chong Peng
//
// Authors: Chong Peng
// Maintainer: Chong Peng and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_singles_casscf_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_singles_casscf_h


#if defined(MPQC_NEW_FEATURES)
#include <chemistry/qc/mbptr12/sr_r12intermediates.h>

namespace sc {

  class CabsSingles {

    public:

    CabsSingles() = default;

    CabsSingles(std::shared_ptr <SingleReference_R12Intermediates<double>> srr12intrmds, bool extra_basis) : singles_r12intrmds_(srr12intrmds), extra_basis_(extra_basis){}

    ~CabsSingles() = default;

    double compute(const std::string &h0);
    void obsolete() { singles_r12intrmds_ = NULL; }
    void print(std::ostream& os = ExEnv::out0()) const;

    const bool extra_basis() {return extra_basis_;}
    const std::shared_ptr <SingleReference_R12Intermediates<double>> r12intermediates() { return singles_r12intrmds_; }


    private:

    template<typename T>
    struct CabsSingles_ {

        typedef TA::Array<T, 4> Array4;
        typedef TA::Array<T, 2> Array2;

        const Array4& Bmatrix;

        CabsSingles_(const Array4& B) : Bmatrix(B){
        }

        /**
         * @param[in] C
         * @param[out] BC
         */
        void operator()(const Array2& C, Array2& BC) {
            BC("x,B'") = Bmatrix("x,B',y,A'") * C("y,A'");
        }
    };

  /// makes a diagonal 2-index preconditioner: pc_x^y = -1/ ( <x|O1|x> - <y|O2|y> )
  template <typename T>
  struct DiagPrecond2 {
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
    DiagPrecond2(const EigenMatrixX& O1_mat,
                  const EigenMatrixX& O2_mat) :
                      O1_mat_(O1_mat), O2_mat_(O2_mat) {
    }
    template <typename Index> T operator()(const Index& i) {
      return 1.0 / (- O1_mat_(i[0], i[0]) + O2_mat_(i[1], i[1]));
    }

    private:
      EigenMatrixX O1_mat_;
      EigenMatrixX O2_mat_;
  };

  std::shared_ptr <SingleReference_R12Intermediates<double>> singles_r12intrmds_;
  bool extra_basis_;

  // compute CABS singles correction using Fock operator as H0
  double CabsSinglesFock();
  // compute CABS singles correction using two-body operators in H0
  double CabsSinglesDyall(const std::string &h0);
  };

}

#endif
#endif // end of header guard CABS_SINGLE

