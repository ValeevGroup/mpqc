//
// transform_tbint.h
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifndef _chemistry_qc_mbptr12_transformtbint_h
#define _chemistry_qc_mbptr12_transformtbint_h

#include <string>
#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>

using namespace std;

namespace sc {

class MOIntsTransformFactory;

  /** TwoBodyMOIntsTransform computes two-body integrals in MO basis
      using parallel integrals-direct AO->MO transformation. */

class TwoBodyMOIntsTransform : virtual public SavableState {
public:
  typedef MOIntsTransformFactory::IntegralCallback IntegralCallback;

private:

  // Construct the integrals accumulator object
  // This function depends on the particulars of the transformation
  virtual void init_acc() = 0;
  // Compute required dynamic memory for a given batch size
  // implementation depends on the particulars of the concrete type
  virtual distsize_t compute_transform_dynamic_memory_(int ni) const = 0;

protected:
  /** By default, integrals smaller than zero_integral are considered zero.
      This constant is only used in checking integrals, not computing them. */
  static const double zero_integral = 1.0e-12;
  /// Predefined enumerated type for the MO spaces
  typedef struct {
    enum {Space1, Space2, Space3, Space4};
  } MOSpaces;

  std::string name_;
  Ref<MOIntsTransformFactory> factory_;

  IntegralCallback callback_;

  Ref<MolecularEnergy> top_mole_;   // Top-level molecular energy to enable checkpointing
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;
  // Integrals accumulator
  Ref<R12IntsAcc> ints_acc_;

  Ref<MOIndexSpace> space1_;
  Ref<MOIndexSpace> space2_;
  Ref<MOIndexSpace> space3_;
  Ref<MOIndexSpace> space4_;

  int num_te_types_;
  size_t memory_;
  bool dynamic_;
  double print_percent_;
  int debug_;
  MOIntsTransformFactory::StoreMethod ints_method_;
  std::string file_prefix_;

  // These variables are never saved but computed every time in case environment
  // has changed or it's a restart  
  size_t mem_static_;
  int batchsize_;
  int npass_;

  /// returns index in range of space1_ where to start the transformation
  unsigned int restart_orbital() const;
  
  // Compute used static memory and batch size
  void init_vars();
  // Re-construct the integrals accumulator object
  void reinit_acc();
  // Allocate distributed memory
  void alloc_mem(const size_t localmem);
  // Deallocate distributed memory
  void dealloc_mem();

  // Compute batchsize given the amount of used static memory and
  // the number of i-orbitals
  int compute_transform_batchsize_(size_t mem_static, int rank_i);

  // Compute the number of ij-pairs per this task
  static int compute_nij(const int rank_i, const int rank_j, const int nproc, const int me);

  /** Generates a report on memory for the transform : user-specified limits, projected and actual use.
      Assumes formatting info from ExEnv::out0().
   */
  void memory_report(std::ostream& os = ExEnv::out0()) const;
  /** Generates a report on MO spaces for the transform.
      Assumes formatting info from ExEnv::out0().
   */
  void mospace_report(std::ostream& os = ExEnv::out0()) const;

  /** Prints out standard header. Call at the beginning of compute().
   */
  void print_header(std::ostream& os = ExEnv::out0()) const;
  /** Prints out standard footer. Call at the end of compute().
   */
  void print_footer(std::ostream& os = ExEnv::out0()) const;

  /** Checks whether this TwoBodyInt is compatible with this TwoBodyMOIntsTransform */
  void check_tbint(const Ref<TwoBodyInt>& tbint) const;

  /// create TwoBodyInt calling callback on integral. Type of params must match callback
  static Ref<TwoBodyInt> create_tbint(const Ref<Integral>& integral,
                                      const IntegralCallback& callback,
                                      const Ref<IntParams>& params);
    
public:

  TwoBodyMOIntsTransform(StateIn&);
  TwoBodyMOIntsTransform(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                         const IntegralCallback& callback,
                         const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                         const Ref<MOIndexSpace>& space3, const Ref<MOIndexSpace>& space4);
  ~TwoBodyMOIntsTransform();

  void save_data_state(StateOut&);

  /// Returns the name of the transform
  std::string name() const {return name_;}
  /// Returns a short label which uniquely identifies the type of transform
  virtual std::string type() const =0;
  /// Returns the MemoryGrp object
  Ref<MemoryGrp> mem() const;
  /// Returns the MessageGrp object
  Ref<MessageGrp> msg() const;
  /** Returns the integrals accumulator object. */
  const Ref<R12IntsAcc>& ints_acc();
  /// Returns MOIndexSpace object 1
  const Ref<MOIndexSpace>& space1() const;
  /// Returns MOIndexSpace object 2
  const Ref<MOIndexSpace>& space2() const;
  /// Returns MOIndexSpace object 3
  const Ref<MOIndexSpace>& space3() const;
  /// Returns MOIndexSpace object 4
  const Ref<MOIndexSpace>& space4() const;

  /// Returns the update print frequency
  double print_percent() const;
  /// Returns the batchsize for each pass of the transformation
  int batchsize() const;
  /// Returns the debug level
  int debug() const;
  /// Returns whether to use dynamic load balancing
  bool dynamic() const;
  /// Returns the number of types of two body integrals computed
  int num_te_types() const;
  /** Returns the number of bytes allocated for each ij-block of integrals of one type
      in MemoryGrp. It's guaranteed to be divisible by sizeof(double).
    */
  virtual const size_t memgrp_blksize() const =0;
  
  /// Specifies the top-level MolecularEnergy object to use for checkpointing
  void set_top_mole(const Ref<MolecularEnergy>& top_mole) { top_mole_ = top_mole; }

  /** Specifies how many integral types computed by TwoBodyInt to be transformed
      Default is 1. */
  void set_num_te_types(const int num_te_types);
  void set_memory(const size_t memory);
  void set_debug(int debug) { debug_ = debug; }
  void set_dynamic(bool dynamic) { dynamic_ = dynamic; }
  void set_print_percent(double print_percent) { print_percent_ = print_percent; }

  /// Computes transformed integrals. param's contents are passed as parameter to the Integral::IntegralCallback
  virtual void compute(const Ref<IntParams>& param) = 0;
  /// Check symmetry of transformed integrals
  virtual void check_int_symm(double threshold = TwoBodyMOIntsTransform::zero_integral) throw (ProgrammingError) =0;
  /// Make the transform obsolete. Next call to compute() will recompute
  virtual void obsolete();

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


