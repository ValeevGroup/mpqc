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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_wfnworld_h
#define _mpqc_src_lib_chemistry_qc_lcao_wfnworld_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/lcao/moints_runtime.h>
#include <chemistry/qc/lcao/fockbuild_runtime.h>

using boost::property_tree::ptree;

namespace sc {

  class XMLWriter;

  /** Class WavefunctionWorld describes the environment of a Wavefunction */
class WavefunctionWorld : virtual public SavableState, virtual public DescribedXMLWritable {

  // change to 0 to use the old set of OrbitalSpace keys
  static const int USE_NEW_ORBITALSPACE_KEYS = 1;

public:

  /// Describes the method of storing transformed MO integrals.
  typedef MOIntsTransform::StoreMethod StoreMethod;

  WavefunctionWorld(StateIn&);
    /** KeyVal constructor uses the following keywords
        <dl>

    <dt><tt>wfn</tt><dd> This specifies the Wavefunction that is in charge of the World ("Czar"). There is no default.

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

    <dt><tt>df</tt><dd> This optional boolean specifies whether to perform density fitting.
    The default is to not perform density fitting, unless <tt>df_basis</tt> is provided.
    @note The default is likely to change before MPQC 3 release.

    <dt><tt>df_basis</tt><dd> This optional GaussianBasisSet object specifies the density-fitting basis
    to use for all density fitting tasks. If <tt>df=true</tt>, the default is to (try to) find the matching
    density fitting basis for the orbital basis.

    <dt><tt>df_kernel</tt><dd> If density fitting is requsted, this keyword will be queried to determine the kernel
    for density-fitting. By default this keyword is not needed as the DensityFittingRuntime and other runtime
    components will try to determine the best density fitting method. If provided, this will provide
    the world-wide choice for the density fitting kernel.
    The only supported values are <tt>coulomb</tt> (this corresponds to the fitting
    the density to reproduce the electric field), <tt>delta</tt> (the overlap), and <tt>exp(X)</tt> (this corresponds to fitting the density to reproduce
    the potential) where <tt>X</tt> is a positive parameter that determines the lengthscale of
    the region in which to fit the potential.

    <dt><tt>df_solver</tt><dd> If df_basis is specified, this keyword will be queried to determine the method by which
    density-fitting will be performed. Valid values are:
    <dl>

      <dt><tt>cholesky_inv</tt><dd> Use Cholesky inverse. Only valid if <tt>df_kernel=coulomb</tt>. This is the cheapest option and the default.
      <dt><tt>cholesky</tt><dd> Use Cholesky linear solver. Only valid if <tt>df_kernel=coulomb</tt>.
      <dt><tt>cholesky_refine</tt><dd> Use Cholesky linear solver + iterative refinement. Only valid if <tt>df_kernel=coulomb</tt>.

      <dt><tt>bunchkaufman_inv</tt><dd> Use Bunch-Kaufman inverse. Valid for any value of <tt>df_kernel</tt>.
      <dt><tt>bunchkaufman</tt><dd> Use Bunch-Kaufman linear solver. Valid for any value of <tt>df_kernel</tt>.
      <dt><tt>bunchkaufman_refine</tt><dd> Use Bunch-Kaufman linear solver + iterative refinement. Valid for any value of <tt>df_kernel</tt>. Use this option
      if you want the maximum numerical precision of the density fitting.
      <dt><tt>householder</tt><dd> Use HouseholderQR decomposition as implemented in Eigen.  Slow for now (doesn't make direct Lapack calls but rather uses Eigen).
      Valid for any value of <tt>df_kernel</tt>.
      <dt><tt>householder_colpiv</tt><dd> Use ColPivHouseholderQR decomposition as implemented in Eigen.  Slow for now (doesn't make direct Lapack calls but rather uses Eigen).
      Valid for any value of <tt>df_kernel</tt>, more numerically stable than plain HouseholderQR
      <dt><tt>householder_fullpiv</tt><dd> Use FullPivHouseholderQR decomposition as implemented in Eigen.  Slow for now (doesn't make direct Lapack calls but rather uses Eigen).
      Valid for any value of <tt>df_kernel</tt>, more numerically stable than plain HouseholderQR or ColPivHouseholderQR

    </dl>

    <dt><tt>df_local_coulomb</tt><dd> If df_basis is specified, this keyword will tell whether strongly local (a.k.a. pair atomic) fitting
    should be used for the coulomb operator

    <dt><tt>df_local_exchange</tt><dd> If df_basis is specified, this keyword will tell whether strongly local (a.k.a. pair atomic) fitting
    should be used for the exchange operator

    <dt><tt>dynamic</tt><dd> This boolean keyword specifies whether dynamic load balancing
    is used by MO integrals transforms. The default is false.

    <dt><tt>ints_precision</tt><dd> This real keyword specifies the precision of MO integrals. Currently, this keyword specifies the precision of AO integrals.
                                    The default is to determine the desired precision heuristically according to
                                    the desired accuracy of the Czar (the current heuristics set the precision to 1e-15 regardless of the Czar accuracy).

        </dl>
    */
  WavefunctionWorld(const Ref<KeyVal>& keyval);
  ~WavefunctionWorld();

  void save_data_state(StateOut&);

  virtual ptree& write_xml(ptree& parent, const XMLWriter& writer);

  /// obsoletes this object
  /// every wavefunction that owns a WavefunctionWorld must call this method when it's obsolete() method is called
  /// @sa Compute::obsolete()
  void obsolete();

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

  Wavefunction* wfn() const { return wfn_; }
  void set_wfn(Wavefunction* w);

  bool df() const { return df_; }
  const Ref<GaussianBasisSet>& basis_df() const { return bs_df_; };

  const Ref<MemoryGrp>& mem() const { return mem_;};
  const Ref<MessageGrp>& msg() const { return msg_;};
  const Ref<ThreadGrp>& thr() const { return thr_;};
  Ref<Integral> integral() const { return wfn()->integral(); }

  bool dynamic() const { return dynamic_; };
  double print_percent() const { return print_percent_; };
  int debug_level() const { return debug_; };
  StoreMethod::type ints_method() const { return ints_method_; };
  const std::string& ints_file() const;
  double ints_precision() const { return ints_precision_; }

  /// Returns the MOIntsTransformFactory object
  const Ref<MOIntsTransformFactory>& tfactory() const { return tfactory_; };
  /// Returns the MOIntsRuntime object
  const Ref<MOIntsRuntime>& moints_runtime() const { return moints_runtime_; };
  //convenient shortcut
  const Ref<TwoBodyFourCenterMOIntsRuntime>& moints_runtime4() const { return moints_runtime_->runtime_4c(); };
  /// Returns the FockBuildRuntime object that can build Fock matrices
  const Ref<FockBuildRuntime>& fockbuild_runtime() const { return fockbuild_runtime_; };


  void print(std::ostream& o) const;

  // call this after obsolete
  void initialize_ao_spaces();

private:

  void xml_data_local(bool do_integrals, ptree& pt, const XMLWriter& writer);
  void xml_data_nonlocal(bool do_integrals, ptree& pt, const XMLWriter& writer);

  /// whose World is it?
  Wavefunction* wfn_;
  bool df_;                       //!< whether to do density fitting
  Ref<GaussianBasisSet> bs_df_;   //!< the density-fitting basis
  std::string df_kernel_;         //!< the density-fitting kernel
  TwoBodyOperSet::type df_kernel_opertype_; //!< operator type of the density-fitting kernel
  Ref<IntParams> df_kernel_params_;     //!< parameters of the density-fitting kernel
  std::string df_solver_;         //!< the density-fitting solver
  bool df_local_coulomb_;         //!< whether or not local density fitting should be employed for the coulomb operator
  bool df_local_exchange_;        //!< whether or not local density fitting should be employed for the exchange operator
  bool exact_diag_J_;             //!< should the exact semidiagonal integrals (ab|ab) be used for coulomb operator?
  bool exact_diag_K_;             //!< should the exact semidiagonal integrals (ab|ab) be used for exchange operator?
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;
  typedef enum {
    DFBasis,
    DFCoefficients,
    DFIntegralsERI,
    ExactIntegralsERI
  } _XMLOutputData;
  std::vector<_XMLOutputData> out_data_;

  bool dynamic_;
  double print_percent_;
  int debug_;
  StoreMethod::type ints_method_;
  std::string ints_file_;
  double ints_precision_;

  /// The transform factory
  Ref<MOIntsTransformFactory> tfactory_;
  /// The MOIntsRuntime object
  Ref<MOIntsRuntime> moints_runtime_;
  /// The FockBuildRuntime object
  Ref<FockBuildRuntime> fockbuild_runtime_;

  void initialize(); // called by constructors

  /// parses kernel key and converts to TwoBodyOper::type and IntParams objects
  static std::pair<TwoBodyOperSet::type, Ref<IntParams> > init_df_kernel(std::string kernel_key);

};

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
