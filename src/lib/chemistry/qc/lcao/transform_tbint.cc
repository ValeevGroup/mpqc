//
// transform_tbint.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <stdexcept>
#include <sstream>
#include <cassert>

#include <util/misc/consumableresources.h>
#include <util/misc/formio.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/lcao/transform_tbint.h>
#include <util/misc/print.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  TwoBodyMOIntsTransform
 -----------*/
static ClassDesc TwoBodyMOIntsTransform_cd(
  typeid(TwoBodyMOIntsTransform),"TwoBodyMOIntsTransform",1,"virtual public SavableState",
  0, 0, 0);

// default values
double TwoBodyMOIntsTransform::zero_integral = 1.0e-12;
int TwoBodyMOIntsTransform::debug_ = 0;
double TwoBodyMOIntsTransform::print_percent_ = 10.0;

TwoBodyMOIntsTransform::TwoBodyMOIntsTransform(const std::string& name, const Ref<MOIntsTransformFactory>& fact,
                                               const Ref<TwoBodyIntDescr>& tbintdescr,
                                               const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                               const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4) :
  name_(name), factory_(fact), tbintdescr_(tbintdescr),
  space1_(space1), space2_(space2), space3_(space3), space4_(space4),
  mem_(factory()->mem()),
  msg_(factory()->msg()),
  thr_(ThreadGrp::get_default_threadgrp()),
  restart_orbital_(0),
  // Default values
  dynamic_(factory()->dynamic()),
  ints_method_(factory()->ints_method()),
  file_prefix_(factory()->file_prefix()),
  max_memory_(ConsumableResources::get_default_instance()->memory()),
  log2_epsilon_(factory_->log2_precision())
{
  // all spaces must be given, even for partial transforms
  MPQC_ASSERT(space1_);
  MPQC_ASSERT(space2_);
  MPQC_ASSERT(space3_);
  MPQC_ASSERT(space4_);
}

TwoBodyMOIntsTransform::TwoBodyMOIntsTransform(StateIn& si) : SavableState(si)
{
  si.get(name_);
  factory_ << SavableState::restore_state(si);
  ints_acc_ << SavableState::restore_state(si);

  space1_ << SavableState::restore_state(si);
  space2_ << SavableState::restore_state(si);
  space3_ << SavableState::restore_state(si);
  space4_ << SavableState::restore_state(si);

  mem_ = MemoryGrp::get_default_memorygrp();
  msg_ = MessageGrp::get_default_messagegrp();
  thr_ = ThreadGrp::get_default_threadgrp();

  int dynamic; si.get(dynamic); dynamic_ = (bool) dynamic;
  int ints_method; si.get(ints_method); ints_method_ = static_cast<MOIntsTransform::StoreMethod::type>(ints_method);
  si.get(file_prefix_);
  si.get(restart_orbital_);
  si.get(log2_epsilon_);
}

TwoBodyMOIntsTransform::~TwoBodyMOIntsTransform()
{
  dealloc_mem();
}

void
TwoBodyMOIntsTransform::save_data_state(StateOut& so)
{
  so.put(name_);
  SavableState::save_state(factory_.pointer(),so);
  SavableState::save_state(ints_acc_.pointer(),so);

  SavableState::save_state(space1_.pointer(),so);
  SavableState::save_state(space2_.pointer(),so);
  SavableState::save_state(space3_.pointer(),so);
  SavableState::save_state(space4_.pointer(),so);

  so.put((int)dynamic_);
  so.put((int)ints_method_);
  so.put(file_prefix_);
  so.put(restart_orbital_);
  so.put(log2_epsilon_);
}

void
TwoBodyMOIntsTransform::set_memgrp(const Ref<MemoryGrp>& new_mem) { mem_ = new_mem; }

const Ref<MemoryGrp>&
TwoBodyMOIntsTransform::mem() const {return mem_; }

const Ref<MessageGrp>&
TwoBodyMOIntsTransform::msg() const {return msg_; }

const Ref<TwoBodyIntDescr>&
TwoBodyMOIntsTransform::intdescr() const { return tbintdescr_; }

const Ref<DistArray4>&
TwoBodyMOIntsTransform::ints_distarray4() {
  init_acc();
  return ints_acc_;
}

void
TwoBodyMOIntsTransform::partially_transformed_ints(const Ref<DistArray4>& ints_acc)
{
  throw FeatureNotImplemented("TwoBodyMOIntsTransform::partially_transformed_ints() -- not implemented for this type of transform",__FILE__,__LINE__);
}

size_t
TwoBodyMOIntsTransform::memory() const {
  return memory_;
}

size_t
TwoBodyMOIntsTransform::peak_memory() const {
  return peak_memory_;
}

const Ref<OrbitalSpace>&
TwoBodyMOIntsTransform::space1() const {return space1_;}

const Ref<OrbitalSpace>&
TwoBodyMOIntsTransform::space2() const {return space2_;}

const Ref<OrbitalSpace>&
TwoBodyMOIntsTransform::space3() const {return space3_;}

const Ref<OrbitalSpace>&
TwoBodyMOIntsTransform::space4() const {return space4_;}

int
TwoBodyMOIntsTransform::batchsize() const {return batchsize_; }

bool
TwoBodyMOIntsTransform::dynamic() const {return dynamic_; }

unsigned int
TwoBodyMOIntsTransform::num_te_types() const { return tbintdescr_->num_sets(); }

unsigned int
TwoBodyMOIntsTransform::restart_orbital() const {
  return restart_orbital_;
}

void
TwoBodyMOIntsTransform::set_log2_epsilon(double prec) {
  if (prec < this->log2_epsilon_)
    this->obsolete();
  log2_epsilon_ = prec;
}

///////////////////////////////////////////////////////
// Compute the batchsize for the transformation
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
///////////////////////////////////////////////////////
int
TwoBodyMOIntsTransform::compute_transform_batchsize_(size_t mem_static, int rank_i)
{
  // Check is have enough for even static objects
  size_t mem_dyn = 0;
  const size_t max_memory = ConsumableResources::get_default_instance()->memory();
  if (max_memory <= mem_static)
    return 0;
  else
    mem_dyn = max_memory - mem_static;

  // Determine if calculation is possible at all (i.e., if ni=1 possible)
  int ni = 1;
  distsize_t maxdyn = compute_transform_dynamic_memory_(ni);
  if (maxdyn > mem_dyn) {
    return 0;
  }

  ni = 2;
  while (ni<=rank_i) {
    maxdyn = compute_transform_dynamic_memory_(ni);
    if (maxdyn >= mem_dyn) {
      ni--;
      break;
    }
    ni++;
  }
  if (ni > rank_i) ni = rank_i;

  return ni;
}


void
TwoBodyMOIntsTransform::init_vars()
{
  const int me = msg_->me();
  const int rank_i = space1_->rank() - restart_orbital_;

  static_memory_ = 0;
  if (me == 0) {
#if 0 // there is not enough information here to figure out how to compute memory requirements -- just add 1 MB
    // mem_static should include storage in OrbitalSpace
    static_memory_ = space1_->memory_in_use() +
                  space2_->memory_in_use() +
                  space3_->memory_in_use() +
                  space4_->memory_in_use(); // scf vector
    int nthreads = thr_->nthread();
    // ... plus the integrals evaluators
    //static_memory_ += nthreads * factory_->integral()->storage_required_grt(space1_->basis(),space2_->basis(),
    //                                                        space3_->basis(),space4_->basis());
#else
    static_memory_ += 0;
#endif
    batchsize_ = compute_transform_batchsize_(static_memory_,rank_i);
  }

  // Send value of ni and mem_static to other nodes
  msg_->bcast(batchsize_);
  double static_memory_double = static_cast<double>(static_memory_);
  msg_->bcast(static_memory_double);
  static_memory_ = static_cast<size_t>(static_memory_double);

  if (batchsize_ == 0)
    throw std::runtime_error("TwoBodyMOIntsTransform::init_vars() -- batch size is 0: more memory or processors are needed");

  npass_ = 0;
  int rest = 0;
  if (batchsize_ == rank_i) {
    npass_ = 1;
    rest = 0;
  }
  else {
    rest = rank_i%batchsize_;
    npass_ = (rank_i - rest)/batchsize_ + 1;
    if (rest == 0) npass_--;
  }

  // At this point I need to figure out how much memory will be used after compute() has been called
  // DistArray4 object will either use none or all of the dynamical memory
  // this will call init_acc() implicitly
  const size_t mem_dyn = distsize_to_size(compute_transform_dynamic_memory_(batchsize_));
  if (!ints_distarray4()->data_persistent()) { // data is held in memory
    memory_ = static_memory_ + mem_dyn;
    peak_memory_ = memory_;
  }
  else { // data is held elsewhere
    memory_ = static_memory_;
    peak_memory_ = memory_ + mem_dyn;
  }
}

void
TwoBodyMOIntsTransform::reinit_acc()
{
  if (ints_acc_)
    ints_acc_ = 0;
  init_acc();
}

void
TwoBodyMOIntsTransform::obsolete()
{
  reinit_acc();
}

void
TwoBodyMOIntsTransform::alloc_mem(const size_t localmem)
{
  if (mem_.null())
    throw std::runtime_error("TwoBodyMOIntsTransform::alloc_mem() -- memory group not initialized");
  mem_->set_localsize(localmem);
  if (debug() >= DefaultPrintThresholds::diagnostics) {
    ExEnv::out0() << indent
                  << "Size of global distributed array:       "
                  << mem_->totalsize()
                  << " Bytes" << endl;
  }
}

void
TwoBodyMOIntsTransform::dealloc_mem()
{
  if (mem_.null())
    throw std::runtime_error("TwoBodyMOIntsTransform::dealloc_mem() -- memory group not initialized");
  mem_->set_localsize(0);
}

int
TwoBodyMOIntsTransform::compute_nij(const int rank_i, const int rank_j, const int nproc, const int me)
{
  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  int index = 0;
  int nij = 0;
  for (int i=0; i<rank_i; i++) {
    for (int j=0; j<rank_j; j++) {
      if (index++ % nproc == 0) nij++;
    }
  }

  return nij;
}

void
TwoBodyMOIntsTransform::memory_report(std::ostream& os) const
{
  const int rank_i_restart = space1_->rank() - restart_orbital_;

  os << indent
     << "Memory available per node:      " << max_memory_ << " Bytes"
     << endl;
  os << indent
     << "Static memory used per node:    " << static_memory_ << " Bytes"
     << endl;
  os << indent
     << "Total memory used per node:     " << peak_memory_ << " Bytes"
     << endl;
  os << indent
     << "Memory required for one pass:   "
     << compute_transform_dynamic_memory_(rank_i_restart)+static_memory_
     << " Bytes"
     << endl;
  os << indent
     << "Minimum memory required:        "
     << compute_transform_dynamic_memory_(1)+static_memory_
     << " Bytes"
     << endl;
  os << indent
     << "Number of passes:               " << (rank_i_restart+batchsize_-1)/batchsize_
     << endl;
  os << indent
     << "Batch size:                     " << batchsize_
     << endl;
}

void
TwoBodyMOIntsTransform::mospace_report(std::ostream& os) const
{
  os << indent << "MO space 1" << endl << incindent;
  space1_->print_summary(os);  os << decindent;
  os << indent << "MO space 2" << endl << incindent;
  space2_->print_summary(os);  os << decindent;
  os << indent << "MO space 3" << endl << incindent;
  space3_->print_summary(os);  os << decindent;
  os << indent << "MO space 4" << endl << incindent;
  space4_->print_summary(os);  os << decindent;
}

void
TwoBodyMOIntsTransform::print_header(std::ostream& os) const
{
  if (debug() >= DefaultPrintThresholds::terse)
    os << indent << "Entered " << name_ << " integrals evaluator (transform type " << type() <<")" << endl;
  os << incindent;

  int nproc = msg_->n();
  if (debug() >= DefaultPrintThresholds::diagnostics)
    os << indent << scprintf("nproc = %i", nproc) << endl;

  if (restart_orbital() && debug() >= DefaultPrintThresholds::diagnostics) {
    os << indent
       << scprintf("Restarting at orbital %d",
                   restart_orbital()) << endl;
  }

  memory_report(os);
  if (dynamic_)
    os << indent << "Using dynamic load balancing." << endl;
  if (debug() >= DefaultPrintThresholds::diagnostics)
    mospace_report(os);
}

void
TwoBodyMOIntsTransform::print_footer(std::ostream& os) const
{
  os << decindent;
  if (debug() >= DefaultPrintThresholds::diagnostics)
    os << indent << "Exited " << name_ << " integrals evaluator (transform type " << type() <<")" << endl;
}

#if 0
void
TwoBodyMOIntsTransform::check_tbint(const Ref<TwoBodyInt>& tbint) const
{
  if (num_te_types_ > tbint->num_tbint_types())
    throw AlgorithmException("TwoBodyMOIntsTransform::check_tbint() -- number of integral types supported by \
current TwoBodyInt is less than\nthe number of types expected by the accumulator",__FILE__,__LINE__);
}
#endif

/////////////////////////////////////////////////////////////////////////////

static ClassDesc TwoBodyThreeCenterMOIntsTransform_cd(
  typeid(TwoBodyThreeCenterMOIntsTransform),"TwoBodyThreeCenterMOIntsTransform",1,
  "virtual public SavableState",
  0, 0, 0);

// default values
int TwoBodyThreeCenterMOIntsTransform::debug_ = 0;
double TwoBodyThreeCenterMOIntsTransform::print_percent_ = 10.0;

TwoBodyThreeCenterMOIntsTransform::~TwoBodyThreeCenterMOIntsTransform() {}

TwoBodyThreeCenterMOIntsTransform::TwoBodyThreeCenterMOIntsTransform(const std::string& name,
  const Ref<MOIntsTransformFactory>& factory,
  const Ref<TwoBodyThreeCenterIntDescr>& tbintdescr,
  const Ref<OrbitalSpace>& space1,
  const Ref<OrbitalSpace>& space2,
  const Ref<OrbitalSpace>& space3) :
    name_(name), factory_(factory), mem_(factory->mem()),
    tbintdescr_(tbintdescr),
    space1_(space1), space2_(space2), space3_(space3),
    ints_method_(factory_->ints_method()),
    file_prefix_(factory_->file_prefix()),
    restart_orbital_(0),
    log2_epsilon_(factory_->log2_precision())
  {
  }

TwoBodyThreeCenterMOIntsTransform::TwoBodyThreeCenterMOIntsTransform(StateIn& si) :
  SavableState(si) {
  si.get(name_);
  factory_ << SavableState::restore_state(si);
  ints_acc_ << SavableState::restore_state(si);
  space1_ << SavableState::restore_state(si);
  space2_ << SavableState::restore_state(si);
  space3_ << SavableState::restore_state(si);

  mem_ = MemoryGrp::get_default_memorygrp();

  int ints_method; si.get(ints_method); ints_method_ = static_cast<MOIntsTransform::StoreMethod::type>(ints_method);
  si.get(file_prefix_);
  si.get(restart_orbital_);
  si.get(log2_epsilon_);
}

void
TwoBodyThreeCenterMOIntsTransform::save_data_state(StateOut& so) {
  so.put(name_);
  SavableState::save_state(factory_.pointer(),so);
  SavableState::save_state(ints_acc_.pointer(),so);
  SavableState::save_state(space1_.pointer(),so);
  SavableState::save_state(space2_.pointer(),so);
  SavableState::save_state(space3_.pointer(),so);

  so.put((int)ints_method_);
  so.put(file_prefix_);
  so.put(restart_orbital_);
  so.put(log2_epsilon_);
}

void
TwoBodyThreeCenterMOIntsTransform::set_memgrp(const Ref<MemoryGrp>& new_mem) { mem_ = new_mem; }

size_t
TwoBodyThreeCenterMOIntsTransform::memory() const {
  return memory_;
}

size_t
TwoBodyThreeCenterMOIntsTransform::peak_memory() const {
  return peak_memory_;
}

unsigned int
TwoBodyThreeCenterMOIntsTransform::num_te_types() const { return tbintdescr_->num_sets(); }

unsigned int
TwoBodyThreeCenterMOIntsTransform::restart_orbital() const {
  return restart_orbital_;
}

void
TwoBodyThreeCenterMOIntsTransform::set_log2_epsilon(double prec) {
  if (prec < this->log2_epsilon_)
    this->obsolete();
  log2_epsilon_ = prec;
}

void
TwoBodyThreeCenterMOIntsTransform::init_vars()
{
  const int me = factory()->msg()->me();
  const int rank_R = space3()->rank() - restart_orbital_;

  static_memory_ = 0;
  if (me == 0) {
#if 0 // there is not enough information here to figure out how to compute memory requirements -- just add 1 MB
    // mem_static should include storage in OrbitalSpace
    static_memory_ = space1()->memory_in_use() +
                  space2()->memory_in_use() +
                  space3()->memory_in_use(); // scf vector
    int nthreads = thr_->nthread();
    // ... plus the integrals evaluators
    //static_memory_ += nthreads * factory_->integral()->storage_required_grt(space1_->basis(),space2_->basis(),
    //                                                        space3_->basis(),space4_->basis());
#else
    static_memory_ += 0;
#endif
  }

  // Send value of ni and mem_static to other nodes
  double static_memory_double = static_cast<double>(static_memory_);
  factory()->msg()->bcast(static_memory_double);
  static_memory_ = static_cast<size_t>(static_memory_double);

  npass_ = 0;
  int rest = 0;
  if (true) {
    npass_ = 1;
    rest = 0;
  }

  // At this point I need to figure out how much memory will be used after compute() has been called
  // DistArray4 object will either use none or all of the dynamical memory
  // this will call init_acc() implicitly
  const size_t mem_dyn = distsize_to_size(compute_transform_dynamic_memory());
  if (!ints_acc()->data_persistent()) { // data is held in memory
    memory_ = static_memory_ + mem_dyn;
    peak_memory_ = memory_;
  }
  else { // data is held elsewhere
    memory_ = static_memory_;
    peak_memory_ = memory_ + mem_dyn;
  }
}

void
TwoBodyThreeCenterMOIntsTransform::reinit_acc()
{
  if (ints_acc_)
    ints_acc_ = 0;
  init_acc();
}

void
TwoBodyThreeCenterMOIntsTransform::obsolete()
{
  reinit_acc();
}

void
TwoBodyThreeCenterMOIntsTransform::alloc_mem(const size_t localmem)
{
  if (mem_.null())
    throw std::runtime_error("TwoBodyThreeCenterMOIntsTransform::alloc_mem() -- memory group not initialized");
  mem_->set_localsize(localmem);
  if (debug() >= DefaultPrintThresholds::diagnostics) {
    ExEnv::out0() << indent
                  << "Size of global distributed array:       "
                  << mem_->totalsize()
                  << " Bytes" << endl;
  }
}

void
TwoBodyThreeCenterMOIntsTransform::dealloc_mem()
{
  if (mem_.null())
    throw std::runtime_error("TwoBodyThreeCenterMOIntsTransform::dealloc_mem() -- memory group not initialized");
  mem_->set_localsize(0);
}

void
TwoBodyThreeCenterMOIntsTransform::memory_report(std::ostream& os) const
{
  os << indent
     << "Memory available per node:      " << max_memory_ << " Bytes"
     << endl;
  os << indent
     << "Static memory used per node:    " << static_memory_ << " Bytes"
     << endl;
  os << indent
     << "Total memory used per node:     " << peak_memory_ << " Bytes"
     << endl;
  this->extra_memory_report(os);
}

void
TwoBodyThreeCenterMOIntsTransform::mospace_report(std::ostream& os) const
{
  os << indent << "MO space 1" << endl << incindent;
  space1_->print_summary(os);  os << decindent;
  os << indent << "MO space 2" << endl << incindent;
  space2_->print_summary(os);  os << decindent;
  os << indent << "MO space 3" << endl << incindent;
  space3_->print_summary(os);  os << decindent;
}

void
TwoBodyThreeCenterMOIntsTransform::print_header(std::ostream& os) const
{
  if (debug() >= DefaultPrintThresholds::terse)
    os << indent << "Entered " << name_ << " integrals evaluator (transform type " << type() <<")" << endl;
  os << incindent;

  int nproc = mem_->n();
  if (debug() >= DefaultPrintThresholds::diagnostics)
    os << indent << scprintf("nproc = %i", nproc) << endl;

  if (restart_orbital() && debug() >= DefaultPrintThresholds::diagnostics) {
    os << indent
       << scprintf("Restarting at orbital %d",
                   restart_orbital()) << endl;
  }

  memory_report(os);
  if (debug() >= DefaultPrintThresholds::diagnostics)
    mospace_report(os);
}

void
TwoBodyThreeCenterMOIntsTransform::print_footer(std::ostream& os) const
{
  os << decindent;
  if (debug() >= DefaultPrintThresholds::diagnostics)
    os << indent << "Exited " << name_ << " integrals evaluator (transform type " << type() <<")" << endl;
}

const Ref<DistArray4>&
TwoBodyThreeCenterMOIntsTransform::ints_acc() {
  init_acc();
  return ints_acc_;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
