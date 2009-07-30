#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>
#include <vector>
#include <ctime>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>

#include <chemistry/qc/ccr12/ccr12.h>

using std::vector;
using namespace std;
using namespace sc;
using namespace sc::LinearR12;

/*--------------------------------
  CCR12
 --------------------------------*/

static ClassDesc CCR12_cd(
  typeid(CCR12),"CCR12",1,"public MBPT2_R12",
  0,create<CCR12>,create<CCR12>);

CCR12::CCR12(StateIn& s): MBPT2_R12(s){
  throw ProgrammingError("sc::CCR12::CCR12(StateIn&) -- constructer not yet implemented",__FILE__,__LINE__);
}


CCR12::CCR12(const Ref<KeyVal>& keyval): MBPT2_R12(keyval){

  keyval_=keyval;


}

void CCR12::common_init(string theory){

  theory_ = theory;
  thrgrp_=ThreadGrp::get_default_threadgrp();
  msggrp_=MessageGrp::get_default_messagegrp();
  mem_=MemoryGrp::get_default_memorygrp();
  timer_=new RegionTimer();

  ExEnv::out0() << endl << indent << "-------- CCR12 calculation --------"<< endl;
  if (sizeof(long)<8)
    ExEnv::out0() << indent << "!!!!!!!!!! \"long int\" less than 8 bytes !!!!!!!!!!" << endl;

  if(theory_=="notheory") throw InputError("CCR12::CCR12 -- no theory specified",__FILE__,__LINE__);
  ExEnv::out0() << endl << indent << "Theory:       " << theory_ << endl << endl;
  perturbative_ = keyval_->stringvalue("perturbative", KeyValValuechar());
  std::transform(perturbative_.begin(), perturbative_.end(), perturbative_.begin(), (int (*)(int))std::toupper);
  ExEnv::out0() << endl << indent << "Perturbative: " << perturbative_ << endl << endl;

  ndiis_=keyval_->intvalue("ndiis",KeyValValueint(2));
  diis_start_=keyval_->intvalue("diis_start",KeyValValueint(0));

  CLSCF* clscfref=dynamic_cast<CLSCF*>(ref().pointer());
  rhf_=(clscfref!=0);

//init_variables(); // inherent in mbpt.h; returns nocc, nsocc, nvir


  //maxiter
  maxiter_=keyval_->intvalue("maxiter",   KeyValValueint(100));
  //cctresh
  ccthresh_=keyval_->doublevalue("ccthresh", KeyValValuedouble(1.0e-9));
  // get the memory sizes
  memorysize_=keyval_->longvalue("memory",   KeyValValueint(200000000));
  ExEnv::out0() << indent << "Memory size per node: " << memorysize_ << endl;
  worksize_=keyval_->longvalue("workmemory",   KeyValValueint(50000000));
  ExEnv::out0() << indent << "Work   size per node: " << worksize_   << endl;
}

CCR12::~CCR12(){
  mem_->sync();
  delete ccr12_info_;
}

void CCR12::compute(){

  r12evalinfo()->initialize();
  r12evalinfo()->set_memory(memorysize_);

  // CCR12_Info will do integral evaluation, before MemoryGrp is used by Tensors
  ccr12_info_=new CCR12_Info(r12evalinfo(),mem_,memorysize_,ref(),nfzcore(),nfzvirt(),
                  molecule()->point_group()->char_table().nirrep(),worksize_,memorysize_,mem_->n(),ndiis_,
                  theory_,perturbative_);

}



/// utilities >>>>>>> form here >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void CCR12::print_iteration_header(string theory){
  if (mem_->me()==0) {
    ExEnv::out0() << endl << indent << theory << " iterations:" << endl;
    ExEnv::out0() << indent << "iter      corr energy        residual RMS        Wall " << endl;
    ExEnv::out0() << indent << "======================================================" << endl;
  }
}

void CCR12::print_iteration_header_short(string theory){
  if (mem_->me()==0) {
    ExEnv::out0() << endl << indent << theory << " iterations:" << endl;
    ExEnv::out0() << indent << "iter      residual RMS        Wall " << endl;
    ExEnv::out0() << indent << "===================================" << endl;
  }
}

void CCR12::print_iteration_footer(){
  if (mem_->me()==0) {
    ExEnv::out0() << indent << "======================================================" << endl;
  }
}

void CCR12::print_iteration_footer_short(){
  if (mem_->me()==0) {
    ExEnv::out0() << indent << "===================================" << endl;
  }
}

void CCR12::print_iteration(int iter,double energy,double rnorm,double iter_start,double iter_end){
  double iter_duration=iter_end-iter_start;
  if (mem_->me()==0) {
    ExEnv::out0() << indent << fixed    << setw(4) << iter
                            << setw(20) << setprecision(13) << energy
                            << setw(20) << setprecision(13) << rnorm
                            << setw(10) << setprecision(2)  << iter_duration << endl;
  }
}

void CCR12::print_iteration_short(int iter,double rnorm,double iter_start,double iter_end){
  double iter_duration=iter_end-iter_start;
  if (mem_->me()==0) {
    ExEnv::out0() << indent << fixed    << setw(4) << iter
                            << setw(20) << setprecision(13) << rnorm
                            << setw(10) << setprecision(2)  << iter_duration << endl;
  }
}

void CCR12::print_timing(double timing, std::string tag){
  if (mem_->me()==0) {
    ExEnv::out0() << indent << fixed    << "Elapsed time [ "<< tag << " ]: "
                            << setw(10) << setprecision(2)  << timing << endl << endl;
  }
}


void CCR12::print_correction(double corr, double base, string theory){
 if (mem_->me()==0) {
  ExEnv::out0() << endl;
  ExEnv::out0() << indent << theory << fixed << setprecision(10) << " correction: " << corr      << endl;
  ExEnv::out0() << indent << theory << fixed << setprecision(10) << " energy    : " << corr+base << endl;
  ExEnv::out0() << endl;
 }
}


void CCR12::print(ostream&o) {
  o << indent << "CCR12:" << endl;
  o << incindent;
  info()->print(o);
  o << decindent;
}
