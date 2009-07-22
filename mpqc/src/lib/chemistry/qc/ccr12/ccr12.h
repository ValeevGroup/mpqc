#ifndef _chemistry_qc_ccr12_ccr12_h
#define _chemistry_qc_ccr12_ccr12_h

#ifdef __GNUC__
#pragma interface
#endif

#include <string>
#include <util/misc/compute.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

//class R12IntEval;
class R12IntEvalInfo;

class CCR12: public MBPT2_R12 {
  protected:
    Ref<ThreadGrp> thrgrp_;
    Ref<MessageGrp> msggrp_;
    Ref<MemoryGrp> mem_;
    Ref<KeyVal> keyval_;
    CCR12_Info* ccr12_info_;
    long worksize_;
    long memorysize_;
    Ref<R12IntEval> r12eval_;
    bool rhf_;
    int maxiter_;
    Ref<RegionTimer> timer_;
    std::string theory_;
    std::string perturbative_;
    double ccthresh_;

    int ndiis_;
    int diis_start_;

  public:
    CCR12(StateIn&);
    CCR12(const Ref<KeyVal>&);
    ~CCR12();
    CCR12_Info* info() const { return ccr12_info_;};


  protected:
    /// always execute this from Derived's compute()
    void compute();
    void common_init(std::string);

    // print nutilities
    void print_iteration_header(std::string);
    void print_iteration_footer();
    void print_iteration(int,double,double,double,double);
    void print_iteration_header_short(std::string);
    void print_iteration_footer_short();
    void print_iteration_short(int,double,double,double);
    void print_correction(double,double,std::string);
    void print(std::ostream&);
    void print_timing(double,std::string);

};

}

#endif
