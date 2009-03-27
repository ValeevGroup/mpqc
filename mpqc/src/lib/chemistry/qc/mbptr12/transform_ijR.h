//
// transform_ijR.h
//
// Copyright (C) 2009 Edward Valeev
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
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_transform_ijR_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_transform_ijR_h

#include <chemistry/qc/basis/intdescr.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>

namespace sc {

  /** TwoBodyMOIntsTransform_ijR computes (ij|R) integrals,
      where R are atomic orbitals, using parallel integral-direct AO->MO transformation.
      */
  class TwoBodyMOIntsTransform_ijR: virtual public SavableState {
    public:
      typedef MOIntsTransformFactory::StorageType StorageType;
      typedef IntegralSetDescr<TwoBodyThreeCenterInt> TwoBodyThreeCenterIntDescr;

      TwoBodyMOIntsTransform_ijR(StateIn&);
      TwoBodyMOIntsTransform_ijR(const std::string& name,
                                 const Ref<MOIntsTransformFactory>& factory,
                                 const Ref<TwoBodyThreeCenterIntDescr>& tbintdescr,
                                 const Ref<OrbitalSpace>& space1,
                                 const Ref<OrbitalSpace>& space2,
                                 const Ref<OrbitalSpace>& space3);
      ~TwoBodyMOIntsTransform_ijR();
      void save_data_state(StateOut&);

      /// Returns the name of the transform
      std::string name() const { return name_; }
      /// Returns a short label which uniquely identifies the type of transform
      std::string type() const {
        return "ijR";
      }

      /// factory who created this
      const Ref<MOIntsTransformFactory>& factory() const { return factory_; }
      /// MemoryGrp object
      const Ref<MemoryGrp>& mem() const { mem_; }
      /// Returns the integral set descriptor
      const Ref<TwoBodyThreeCenterIntDescr>& intdescr() const { return tbintdescr_; }
      /** Returns the integrals accumulator object. */
      const Ref<R12IntsAcc>& ints_acc();
      /// Returns OrbitalSpace object 1
      const Ref<OrbitalSpace>& space1() const { return space1_; }
      /// Returns OrbitalSpace object 2
      const Ref<OrbitalSpace>& space2() const { return space2_; }
      /// Returns OrbitalSpace object 3
      const Ref<OrbitalSpace>& space3() const { return space3_; }

      /// Returns amount of memory used by this object after compute() has been called
      size_t memory() const;
      /// Returns the maximum amount of memory that will be used by this object
      size_t peak_memory() const;

      /// Returns the update print frequency
      double print_percent() const;
      /// Returns the batchsize for each pass of the transformation
      int batchsize() const;
      /// Returns the debug level
      int debug() const;
      /// Returns the number of types of two body integrals computed
      unsigned int num_te_types() const;

      /// Computes transformed integrals
      void compute();
      /// Make the transform obsolete. Next call to compute() will recompute
      void obsolete();

    private:

      std::string name_;
      Ref<MOIntsTransformFactory> factory_;
      Ref<MemoryGrp> mem_;
      Ref<TwoBodyThreeCenterIntDescr> tbintdescr_;
      // Integrals accumulator
      Ref<R12IntsAcc> ints_acc_;

      Ref<OrbitalSpace> space1_;
      Ref<OrbitalSpace> space2_;
      Ref<OrbitalSpace> space3_;

      int restart_orbital_; // when restarting, this recalls where to start transform
      size_t peak_memory_;  // actual maximum memory (per process) used by this transform during its lifetime
      size_t memory_;       // memory (per process) used by this transform after compute() has been called
      double print_percent_;
      int debug_;
      MOIntsTransformFactory::StoreMethod::type ints_method_;
      std::string file_prefix_;

      // These variables computed every time in case environment has changed or it's a restart
      size_t max_memory_;      // max memory given to this object
      size_t static_memory_;   // memory used to hold persistent quantities
      int batchsize_;
      int npass_;

      /// returns index in range of space1_ where to start the transformation
      unsigned int restart_orbital() const;

      // Compute used static memory and batch size
      void init_vars();
      // Construct the integrals accumulator object
      void init_acc();
      // Re-construct the integrals accumulator object
      void reinit_acc();
      // Allocate distributed memory
      void alloc_mem(const size_t localmem);
      // Deallocate distributed memory
      void dealloc_mem();
      // Reset mem_ to point to new_mem
      void set_memgrp(const Ref<MemoryGrp>& new_mem);

      // Compute required dynamic memory for a given batch size
      // implementation depends on the particulars of the concrete type
      distsize_t compute_transform_dynamic_memory(int nR) const;

      // Compute batchsize given the amount of used static memory and
      // the number of i-orbitals
      int compute_transform_batchsize(size_t mem_static, int rank_R);

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

  };

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
