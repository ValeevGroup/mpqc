//
// wfnworld.h
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

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_wfnworld_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_wfnworld_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/mbptr12/moints_runtime.h>
#include <chemistry/qc/mbptr12/fockbuild_runtime.h>

namespace sc {

  /** Class WavefunctionWorld describes the environment of a Wavefunction */
class WavefunctionWorld : virtual public SavableState {

  // change to 0 to use the old set of OrbitalSpace keys
  static const int USE_NEW_ORBITALSPACE_KEYS = 1;

public:

  /// Describes the method of storing transformed MO integrals.
  typedef MOIntsTransform::StoreMethod StoreMethod;

  WavefunctionWorld(StateIn&);
    /** KeyVal constructor uses the following keywords
        <dl>

    <dt><tt>store_ints</tt><dd> This specifies how to store transformed MO integrals.
    Valid values are:

    <dl>

      <dt><tt>posix</tt><dd> Store integrals in a binary file on task 0's node using POSIX I/O.
      This method does not allow all steps to be parallelized but it is most likely to work in all environments.

      <dt><tt>mpi</tt><dd> Store integrals in a binary file using MPI-I/O. This method allows
      parallelization of all steps, but requires MPI-I/O capability (including MPI-I/O capable file system;
      see keyword <tt>ints_file</tt>)

      <dt><tt>mem</tt><dd> Store integrals in memory. Can only be used with single-pass
      transformations for MP2-R12/A and MP2-R12/A' methods. This choice is the most efficient, but
      requires significant amount memory. It is probably only feasible when <tt>stdapprox = A'</tt>.

      <dt><tt>mem-posix</tt><dd> The program will choose between <tt>mem</tt> and <tt>posix</tt> automatically.

      <dt><tt>mem-mpi</tt><dd> The program will choose between <tt>mem</tt> and <tt>mpi</tt> automatically.

    </dl>

    The default is <tt>posix</tt>.

    <dt><tt>ints_file</tt><dd> This specifies the prefix for the transformed
    MO integrals file if <tt>ints</tt> is set to <tt>posix</tt> or <tt>mpi</tt>.
    If the prefix ends in '/' (slash character) then <i>basename</i><tt>.moints</tt>
    is appended to it where <i>basename</i> is the basename as defined in SCFormIO.
    The default value for the prefix is "./".
    If MPI-I/O is used then it is user's responsibility to ensure
    that the file resides on a file system that supports MPI-I/O.

    <dt><tt>dynamic</tt><dd> This boolean keyword specifies whether dynamic load balancing
    is used by MO integrals transforms. The default is false.


        </dl>
    */
  WavefunctionWorld(const Ref<KeyVal>& keyval,
                    Wavefunction* wfn);
  ~WavefunctionWorld();

  void save_data_state(StateOut&);

  /** Sets whether to use dynamic load balancing in parallel MO transformations. */
  void dynamic(bool dynamic) { dynamic_ = dynamic; };
  /// Sets how frequently updates of progress are printed out. Default is 10%
  void print_percent(double print_percent) { print_percent_ = print_percent; };
  /// Set debug level. Default is 0.
  void debug_level(int debug) { debug_ = debug; };
  /** Sets the method of storing transformed MO integrals. Default depends on
      how the object was constructed. */
  void ints_method(const StoreMethod::type method) { ints_method_ = method; };
  /** Sets name of the file used to store transformed integrals.
      Default depends on how the object was constructed. */
  void ints_file(const std::string& filename) { ints_file_ = filename; };
  /** Sets the amount of memory to use for the calculation. Default is
      determined by DEFAULT_SC_MEMORY. */
  void memory(const size_t nbytes);

  Wavefunction* wfn() const { return wfn_; }

  const Ref<GaussianBasisSet>& basis_df() const { return bs_df_; };

  const Ref<MemoryGrp>& mem() const { return mem_;};
  const Ref<MessageGrp>& msg() const { return msg_;};
  const Ref<ThreadGrp>& thr() const { return thr_;};
  Ref<Integral> integral() const { return wfn()->integral(); }

  bool dynamic() const { return dynamic_; };
  double print_percent() const { return print_percent_; };
  int debug_level() const { return debug_; };
  const StoreMethod::type ints_method() const { return ints_method_; };
  const std::string& ints_file() const;
  const size_t memory() const { return memory_; };

  /// Returns the MOIntsTransformFactory object
  const Ref<MOIntsTransformFactory>& tfactory() const { return tfactory_; };
  /// Returns the MOIntsRuntime object
  const Ref<MOIntsRuntime>& moints_runtime() const { return moints_runtime_; };
  //convenient shortcut
  const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_runtime4() const { return moints_runtime_->runtime_4c(); };
  /// Returns the MOIntsRuntime object
  const Ref<FockBuildRuntime>& fockbuild_runtime() const { return fockbuild_runtime_; };

  void print(std::ostream& o) const;

private:

  /// whose World is it?
  Wavefunction* wfn_;
  Ref<GaussianBasisSet> bs_df_;   //!< the density-fitting basis
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;

  size_t memory_;
  bool dynamic_;
  double print_percent_;
  int debug_;
  StoreMethod::type ints_method_;
  std::string ints_file_;

  /// The transform factory
  Ref<MOIntsTransformFactory> tfactory_;
  /// The MOIntsRuntime object
  Ref<MOIntsRuntime> moints_runtime_;
  /// The FockBuildRuntime object
  Ref<FockBuildRuntime> fockbuild_runtime_;

  void initialize(); // called by constructors
};

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
