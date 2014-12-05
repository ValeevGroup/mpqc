//
// cabs_single.h
//
#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_singles_casscf_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_singles_casscf_h


#if defined(HAVE_MPQC3_RUNTIME)
#include <chemistry/qc/mbptr12/sr_r12intermediates.h>


namespace sc {

  class CabsSingles {

    public:
    CabsSingles(std::shared_ptr <SingleReference_R12Intermediates<double>> srr12intrmds);

    double compute(const std::string &h0);

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

    std::shared_ptr <SingleReference_R12Intermediates<double>> srr12intrmds_;

    // compute CABS singles correction using Fock operator as H0
    double CabsSinglesFock();
    // compute CABS singles correction using two-body operators in H0
    double CabsSinglesDyall(const std::string &h0);
  };

}

#endif
#endif // end of header guard CABS_SINGLE

