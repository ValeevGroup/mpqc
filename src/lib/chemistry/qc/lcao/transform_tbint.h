//
// transform_tbint.h
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

#ifndef _chemistry_qc_lcao_transformtbint_h
#define _chemistry_qc_lcao_transformtbint_h

#include <string>
#include <util/ref/ref.h>
#include <util/misc/scexception.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/intdescr.h>
#include <chemistry/qc/basis/distshpair.h>
#include <math/distarray4/distarray4.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/lcao/transform_factory.h>

namespace sc {

  /** TwoBodyMOIntsTransform computes two-body integrals in MO basis
      using parallel integrals-direct AO->MO transformation.

      The target MO integrals are put into an DistArray4 object.

      */
class TwoBodyMOIntsTransform : virtual public SavableState {
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
  static double zero_integral;
  /// Predefined enumerated type for the MO spaces
  typedef struct {
    enum {Space1, Space2, Space3, Space4};
  } MOSpaces;

  std::string name_;
  Ref<MOIntsTransformFactory> factory_;

  Ref<MolecularEnergy> top_mole_;   // Top-level molecular energy to enable checkpointing
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;    // Initially taken from Factory, but may decide to create another MemoryGrp object
                          // e.g. MemoryGrpRegion object created from Factory's MemoryGrp object. The choice of MemoryGrp
                          // is finialized in the derived class constructor (\sa init_acc())
  Ref<ThreadGrp> thr_;
  Ref<TwoBodyIntDescr> tbintdescr_;
  // Integrals accumulator
  Ref<DistArray4> ints_acc_;

  Ref<OrbitalSpace> space1_;
  Ref<OrbitalSpace> space2_;
  Ref<OrbitalSpace> space3_;
  Ref<OrbitalSpace> space4_;

  int restart_orbital_; // when restarting, this recalls where to start transform
  size_t peak_memory_;  // actual maximum memory (per process) used by this transform during its lifetime
  size_t memory_;       // memory (per process) used by this transform after compute() has been called
  bool dynamic_;
  DistShellPair::SharedData spdata_;
  MOIntsTransform::StoreMethod::type ints_method_;
  std::string file_prefix_;
  double log2_epsilon_;

  // These variables computed every time in case environment has changed or it's a restart
  size_t max_memory_;      // max memory given to this object
  size_t static_memory_;   // memory used to hold persistent quantities
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
  // Reset mem_ to point to new_mem
  void set_memgrp(const Ref<MemoryGrp>& new_mem);

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

#if 0
  /** Checks whether this TwoBodyInt is compatible with this TwoBodyMOIntsTransform */
  void check_tbint(const Ref<TwoBodyInt>& tbint) const;
#endif

public:

  TwoBodyMOIntsTransform(StateIn&);
  TwoBodyMOIntsTransform(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                         const Ref<TwoBodyIntDescr>& tbintdescr,
                         const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                         const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4);
  virtual ~TwoBodyMOIntsTransform();

  void save_data_state(StateOut&);

  /// factory who created this
  const Ref<MOIntsTransformFactory>& factory() const { return factory_; }
  /// Returns the name of the transform
  std::string name() const {return name_;}
  /// Returns a short label which uniquely identifies the type of transform
  virtual std::string type() const =0;
  /// Returns the MemoryGrp object
  const Ref<MemoryGrp>& mem() const;
  /// Returns the MessageGrp object
  const Ref<MessageGrp>& msg() const;
  /// Returns the integral set descriptor
  const Ref<TwoBodyIntDescr>& intdescr() const;
  /** Returns the DistArray4 object that holds the integrals. */
  const Ref<DistArray4>& ints_distarray4();
  /// Returns OrbitalSpace object 1
  const Ref<OrbitalSpace>& space1() const;
  /// Returns OrbitalSpace object 2
  const Ref<OrbitalSpace>& space2() const;
  /// Returns OrbitalSpace object 3
  const Ref<OrbitalSpace>& space3() const;
  /// Returns OrbitalSpace object 4
  const Ref<OrbitalSpace>& space4() const;
  /**
   *
   * @param \f$ \log_2(\epsilon) \f$, where \f$ \epsilon \f$ is the absolute numerical precision of the integrals
   *        computed by this object.
   */
  double log2_epsilon() const { return log2_epsilon_; }
  /// \sa log2_epsilon()
  void set_log2_epsilon(double prec);

  /// Supplies the partially transformed integrals.
  virtual void partially_transformed_ints(const Ref<DistArray4>&);


  /// Returns amount of memory used by this object after compute() has been called
  size_t memory() const;
  /// Returns the maximum amount of memory that will be used by this object
  size_t peak_memory() const;

  /// Returns the batchsize for each pass of the transformation
  int batchsize() const;
  /// Returns whether to use dynamic load balancing
  bool dynamic() const;
  /// Returns the number of types of two body integrals computed
  unsigned int num_te_types() const;
  /** Returns the number of bytes allocated for each ij-block of integrals of one type
      in MemoryGrp. It's guaranteed to be divisible by sizeof(double).
    */
  virtual size_t memgrp_blksize() const =0;

  /// Specifies the top-level MolecularEnergy object to use for checkpointing
  void set_top_mole(const Ref<MolecularEnergy>& top_mole) { top_mole_ = top_mole; }
  void set_dynamic(bool dynamic) { dynamic_ = dynamic; }

  /// Computes transformed integrals
  virtual void compute() = 0;
  /// Check symmetry of transformed integrals
  virtual void check_int_symm(double threshold = TwoBodyMOIntsTransform::zero_integral) =0;
  /// Make the transform obsolete. Next call to compute() will recompute
  virtual void obsolete();

  /** Returns a that data that must be shared between all DistShellPair
   * objects. */
  DistShellPair::SharedData *shell_pair_data() { return &spdata_; }

public:
  static void set_print_percent(double pp) { print_percent_ = pp; }
  static double print_percent() { return print_percent_; }
  static void set_debug(int d) { debug_ = d; }
  static int debug() { return debug_; }

private:
  static double print_percent_;
  static int debug_;

};


/** TwoBodyThreeCenterMOIntsTransform computes (xy|z) integrals,
    using parallel integral-direct AO->MO transformation.

    The target MO integrals are put into an DistArray4 object stored as (0 x|y z), where 0 is the
    dummy index.

    */
class TwoBodyThreeCenterMOIntsTransform: virtual public SavableState {
  public:
    TwoBodyThreeCenterMOIntsTransform(StateIn&);
    TwoBodyThreeCenterMOIntsTransform(const std::string& name,
                               const Ref<MOIntsTransformFactory>& factory,
                               const Ref<TwoBodyThreeCenterIntDescr>& tbintdescr,
                               const Ref<OrbitalSpace>& space1,
                               const Ref<OrbitalSpace>& space2,
                               const Ref<OrbitalSpace>& space3);
    ~TwoBodyThreeCenterMOIntsTransform();
    void save_data_state(StateOut&);

    /// Returns the name of the transform
    std::string name() const { return name_; }
    /// Returns a short label which uniquely identifies the type of transform
    virtual std::string type() const =0;

    /// factory who created this
    const Ref<MOIntsTransformFactory>& factory() const { return factory_; }
    /// MemoryGrp object
    const Ref<MemoryGrp>& mem() const { return mem_; }
    /// Returns the integral set descriptor
    const Ref<TwoBodyThreeCenterIntDescr>& intdescr() const { return tbintdescr_; }
    /** Returns the integrals accumulator object. */
    const Ref<DistArray4>& ints_acc();
    /// Returns OrbitalSpace object 1
    const Ref<OrbitalSpace>& space1() const { return space1_; }
    /// Returns OrbitalSpace object 2
    const Ref<OrbitalSpace>& space2() const { return space2_; }
    /// Returns OrbitalSpace object 3
    const Ref<OrbitalSpace>& space3() const { return space3_; }
    /**
     *
     * @param \f$ \log_2(\epsilon) \f$, where \f$ \epsilon \f$ is the absolute numerical precision of the integrals
     *        computed by this object.
     */
    double log2_epsilon() const { return log2_epsilon_; }
    /// \sa log2_epsilon()
    void set_log2_epsilon(double prec);

    /// Returns amount of memory used by this object after compute() has been called
    size_t memory() const;
    /// Returns the maximum amount of memory that will be used by this object
    size_t peak_memory() const;
    /// Returns the number of types of two body integrals computed
    unsigned int num_te_types() const;

    /// Computes transformed integrals
    virtual void compute() =0;
    /// Make the transform obsolete. Next call to compute() will recompute
    void obsolete();

  protected:

    std::string name_;
    Ref<MOIntsTransformFactory> factory_;
    Ref<MemoryGrp> mem_;
    Ref<TwoBodyThreeCenterIntDescr> tbintdescr_;
    // Integrals accumulator
    Ref<DistArray4> ints_acc_;

    Ref<OrbitalSpace> space1_;
    Ref<OrbitalSpace> space2_;
    Ref<OrbitalSpace> space3_;

    int restart_orbital_; // when restarting, this recalls where to start transform
    size_t peak_memory_;  // actual maximum memory (per process) used by this transform during its lifetime
    size_t memory_;       // memory (per process) used by this transform after compute() has been called
    MOIntsTransform::StoreMethod::type ints_method_;
    std::string file_prefix_;
    double log2_epsilon_;

    // These variables computed every time in case environment has changed or it's a restart
    size_t max_memory_;      // max memory given to this object
    size_t static_memory_;   // memory used to hold persistent quantities
    int batchsize_;
    int npass_;

    /// returns index in range of space1_ where to start the transformation
    unsigned int restart_orbital() const;

    // Compute used static memory and batch size
    virtual void init_vars();
    // Construct the integrals accumulator object
    virtual void init_acc() =0;
    // Re-construct the integrals accumulator object
    void reinit_acc();
    // Allocate distributed memory
    void alloc_mem(const size_t localmem);
    // Deallocate distributed memory
    void dealloc_mem();
    // Reset mem_ to point to new_mem
    void set_memgrp(const Ref<MemoryGrp>& new_mem);

    /** Compute required dynamic memory for the given batchsize.
        batchsize=-1 means use the default batchsize = space3()->rank()
     */
    virtual distsize_t compute_transform_dynamic_memory(int batchsize = -1) const =0;

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

  public:
    static void set_print_percent(double pp) { print_percent_ = pp; }
    static double print_percent() { return print_percent_; }
    static void set_debug(int d) { debug_ = d; }
    static int debug() { return debug_; }

  private:
    static double print_percent_;
    static int debug_;

    /** Extra memory report. */
    virtual void extra_memory_report(std::ostream& os = ExEnv::out0()) const =0;

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


