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
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>

using namespace std;

namespace sc {

class MOIntsTransformFactory;

  /** TwoBodyMOIntsTransform computes two-body integrals in MO basis
      using parallel integrals-direct AO->MO transformation. */

class TwoBodyMOIntsTransform : virtual public SavableState {

protected:
  std::string name_;
  Ref<MOIntsTransformFactory> factory_;

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

  // Other cases will be handled later
  static const int num_te_types_ = 3;

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
  
  // Compute used static memory and batch size
  void init_vars();
  // Construct the integrals accumulator object
  virtual void init_acc() = 0;
  // Re-construct the integrals accumulator object
  void reinit_acc();

  // Compute batchsize given the amount of used static memory and
  // the number of i-orbitals
  int compute_transform_batchsize_(size_t mem_static, int rank_i);
  
  // Compute required dynamic memory for a given batch size
  // implementation depends on the particulars of the concrete type
  virtual distsize_t compute_transform_dynamic_memory_(int ni) const = 0;

public:

  TwoBodyMOIntsTransform(StateIn&);
  TwoBodyMOIntsTransform(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                         const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                         const Ref<MOIndexSpace>& space3, const Ref<MOIndexSpace>& space4);
  ~TwoBodyMOIntsTransform();

  void save_data_state(StateOut&);

  /// Specifies the top-level MolecularEnergy object to use for checkpointing
  void set_top_mole(const Ref<MolecularEnergy>& top_mole) { top_mole_ = top_mole; }

  void set_debug(int debug) { debug_ = debug; }
  void set_dynamic(bool dynamic) { dynamic_ = dynamic; }
  void set_print_percent(double print_percent) { print_percent_ = print_percent; }

  /// Returns MOIndexSpace object 1
  Ref<MOIndexSpace> space1() const;
  /// Returns MOIndexSpace object 2
  Ref<MOIndexSpace> space2() const;
  /// Returns MOIndexSpace object 3
  Ref<MOIndexSpace> space3() const;
  /// Returns MOIndexSpace object 4
  Ref<MOIndexSpace> space4() const;

  /// Computes transformed integrals
  virtual void compute() = 0;
  
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


