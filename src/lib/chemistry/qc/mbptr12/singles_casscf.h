//
// cabs_single.h
//
#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_singles_casscf_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_singles_casscf_h


#if defined(MPQC_NEW_FEATURES)
#include <chemistry/qc/mbptr12/sr_r12intermediates.h>


namespace sc {

  class CABS_Single {

    public:
    CABS_Single(std::shared_ptr <SingleReference_R12Intermediates<double>> srr12intrmds);

    double compute(const std::string &h0);

    private:
    std::shared_ptr <SingleReference_R12Intermediates<double>> srr12intrmds_;

    // compute CABS singles correction using Fock operator as H0
    double cabs_singles_Fock();
    // compute CABS singles correction using two-body operators in H0
    double cabs_singles_Dyall(const std::string &h0);
  };

}

#endif
#endif // end of header guard CABS_SINGLE

