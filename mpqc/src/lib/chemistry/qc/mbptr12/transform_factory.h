//
// transform_factory.h
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

#ifndef _chemistry_qc_mbptr12_transformfactory_h
#define _chemistry_qc_mbptr12_transformfactory_h

#include <string>
#include <util/ref/ref.h>
#include <util/group/memory.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/linearr12.h>

using namespace std;

namespace sc {

class TwoBodyMOIntsTransform;
//class TwoBodyIntDescr;

  /** MOIntsTransformFactory is a factory that produces MOIntsTransform objects. */

class MOIntsTransformFactory : virtual public SavableState {

public:

  /// Describes the method of storing transformed MO integrals.
  struct StoreMethod {
    enum type { mem_posix = 0, posix = 1, mem_mpi = 2, mpi = 3, mem_only = 4 };
  };
  /// How integrals are stored. Type_13 means (ix|jy) integrals are stored as (ij|xy)
  enum StorageType {StorageType_First=0, StorageType_Last=1,
                    StorageType_12=0, StorageType_13=1};

private:
  
  /// The default integral descriptor will compute ERIs
  typedef TwoBodyIntDescrERI DefaultTwoBodyIntDescr;
  
  Ref<MolecularEnergy> top_mole_;   // Top-level molecular energy to enable checkpointing

  Ref<Integral> integral_;
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;
  
  Ref<TwoBodyIntDescr> tbintdescr_;
  
  Ref<MOIndexSpace> space1_;
  Ref<MOIndexSpace> space2_;
  Ref<MOIndexSpace> space3_;
  Ref<MOIndexSpace> space4_;

  size_t memory_;
  bool dynamic_;
  double print_percent_;
  int debug_;
  StoreMethod::type ints_method_;
  std::string file_prefix_;

public:

  MOIntsTransformFactory(StateIn&);
  MOIntsTransformFactory(const Ref<Integral>& integral,
                         const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2 = 0,
                         const Ref<MOIndexSpace>& space3 = 0, const Ref<MOIndexSpace>& space4 = 0);
  ~MOIntsTransformFactory();

  void save_data_state(StateOut&);

  /// Sets the orbital spaces
  void set_spaces(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2 = 0,
                  const Ref<MOIndexSpace>& space3 = 0, const Ref<MOIndexSpace>& space4 = 0);

  /// Specifies the top-level MolecularEnergy object to use for checkpointing
  void set_top_mole(const Ref<MolecularEnergy>& top_mole) { top_mole_ = top_mole; }
  /// Changes the default TwoBodyIntDescr used to produce integrals
  void tbintdescr(const Ref<TwoBodyIntDescr>& descr) { tbintdescr_ = descr; }
  /// Sets the method of storing transformed MO integrals. Default method is mem_posix.
  void set_ints_method(const StoreMethod::type method) { ints_method_ = method; }
  /// Sets the name of the file to hold the integrals.
  void set_file_prefix(const std::string& prefix) { file_prefix_ = prefix; }
  void set_debug(int debug) { debug_ = debug; }
  void set_dynamic(bool dynamic) { dynamic_ = dynamic; }
  void set_print_percent(double print_percent) { print_percent_ = print_percent; }
  void set_memory(size_t nbytes) { memory_ = nbytes; }

  /// Returns the Integral factory
  Ref<Integral> integral() const { return integral_; };
  /// Returns the default TwoBodyIntDescr used to produce integrals
  Ref<TwoBodyIntDescr> tbintdescr() const { return tbintdescr_; }
  /// Returns the method of storing transformed MO integrals.
  const StoreMethod::type ints_method() const { return ints_method_; }
  /// Sets the name of the file to hold the integrals.
  const std::string file_prefix() const { return file_prefix_; }
  const int debug() const { return debug_; }
  const bool dynamic() const { return dynamic_; }
  const double print_percent() const { return print_percent_; }
  const size_t memory() const { return memory_; }

  /// Returns MOIndexSpace object 1
  Ref<MOIndexSpace> space1() const;
  /// Returns MOIndexSpace object 2
  Ref<MOIndexSpace> space2() const;
  /// Returns MOIndexSpace object 3
  Ref<MOIndexSpace> space3() const;
  /// Returns MOIndexSpace object 4
  Ref<MOIndexSpace> space4() const;

  /** Creates an TwoBodyMOIntsTransform object that will compute (pq|rs) integrals
      stored in qs blocks for each pr */
  Ref<TwoBodyMOIntsTransform>
  twobody_transform_13(const std::string& id, const Ref<TwoBodyIntDescr>& descr = 0);

  /** Creates an TwoBodyMOIntsTransform object that will compute (pq|rs) integrals
    stored in rs blocks for each pq */
  Ref<TwoBodyMOIntsTransform>
  twobody_transform_12(const std::string& id, const Ref<TwoBodyIntDescr>& descr = 0);
  
  /** Creates an TwoBodyMOIntsTransform object that will compute (pq|rs) integrals
    stored according to storage */
  Ref<TwoBodyMOIntsTransform>
  twobody_transform(StorageType storage, const std::string& id,
                    const Ref<TwoBodyIntDescr>& descr = 0);
};

}

#include <chemistry/qc/mbptr12/transform_tbint.h>

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


