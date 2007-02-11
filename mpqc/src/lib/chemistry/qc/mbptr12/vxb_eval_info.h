//
// vxb_eval_info.h
//
// Copyright (C) 2003 Edward Valeev
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

#include <string>
#include <util/misc/string.h>
#include <util/ref/ref.h>
#include <math/scmat/abstract.h>
#include <util/group/memory.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/ansatz.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#include <chemistry/qc/mbptr12/singlerefinfo.h>

#ifndef _chemistry_qc_mbptr12_vxbevalinfo_h
#define _chemistry_qc_mbptr12_vxbevalinfo_h

namespace sc {

class MBPT2_R12;

  /** Class R12IntEvalInfo contains information necessary for R12 intermediate
      evaluators */

class R12IntEvalInfo : virtual public SavableState {

public:

  /// Describes the method of storing transformed MO integrals. See MBPT2_R12.
  typedef MOIntsTransformFactory::StoreMethod StoreMethod;

  /// Maintains virtual orbitals and RI space info if VBS != OBS
  typedef struct {
    //Ref<MOIndexSpace> vbs_sb_;
    //Ref<MOIndexSpace> vbs_;
    Ref<MOIndexSpace> vir_sb_;
    Ref<MOIndexSpace> vir_;
    Ref<MOIndexSpace> vir_act_;
    // RI space
    Ref<MOIndexSpace> ri_;
    // "constructor" that uses SingleRefInfo object
    void init(const Ref<SingleRefInfo>& refinfo, const SpinCase1& spincase);
    // "constructor" that uses SingleRefInfo object
    void init(const Ref<SingleRefInfo>& refinfo, const SpinCase1& spincase, const Ref<MOIndexSpace>& vbs);
  } SpinSpaces;

private:

  Wavefunction* wfn_;     // Wave function that owns this
  Ref<GaussianBasisSet> bs_aux_;
  Ref<GaussianBasisSet> bs_vir_;
  Ref<GaussianBasisSet> bs_ri_;
  Ref<SCMatrixKit> matrixkit_;
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;

  size_t memory_;
  bool dynamic_;
  double print_percent_;
  int debug_;

  Ref<LinearR12::CorrelationFactor> corrfactor_;
  LinearR12::StandardApproximation stdapprox_;
  Ref<LinearR12Ansatz> ansatz_;
  LinearR12::ABSMethod abs_method_;
  bool gbc_;
  bool ebc_;
  bool spinadapted_;
  bool include_mp1_;
  StoreMethod::type ints_method_;
  std::string ints_file_;
  unsigned int maxnabs_;
  /// should I follow Klopper-Samson approach in the intermediates formulation for the EBC-free method?
  bool ks_ebcfree_;
  /// omit P from intermediate B under standard approximation B? (matters only if OBS != ABS)
  bool omit_P_;
  bool safety_check_;
  LinearR12::PositiveDefiniteB posdef_B_;

  int nlindep_aux_;
  int nlindep_vir_;
  int nlindep_ri_;
  
  /// This space depends only on the orbital basis set
  Ref<MOIndexSpace> abs_space_;  // ABS space
  Ref<MOIndexSpace> ribs_space_; // RIBS basis  
  SpinSpaces vir_spaces_[NSpinCases1];
  Ref<MOIndexSpace> vir_act_;
  Ref<MOIndexSpace> vir_;
  Ref<MOIndexSpace> vir_sb_;
  /// Initializes all spaces that relate to the reference determinant
  Ref<SingleRefInfo> refinfo_;
  
  /// The transform factory
  Ref<MOIntsTransformFactory> tfactory_;

  // construct the RI basis based on abs_method
  void construct_ri_basis_(bool safe);
  void construct_ri_basis_ks_(bool safe);
  void construct_ri_basis_ksplus_(bool safe);
  void construct_ri_basis_ev_(bool safe);
  void construct_ri_basis_evplus_(bool safe);
  // Uses ri_basis to construct a basis that spans the orthogonal complement to the OBS
  void construct_ortho_comp_svd_();
  // Returns true if ABS spans OBS
  bool abs_spans_obs_();
  // Construct orthog_aux_
  void construct_orthog_aux_();
  // Construct orthog_vir_
  void construct_orthog_vir_();
  // Construct orthog_ri_
  void construct_orthog_ri_();
  // Throw ProgrammingError if reference is spin_polarized()
  void throw_if_spin_polarized() const;

public:
  R12IntEvalInfo(StateIn&);
  /// Constructs an R12IntEvalInfo object using data from the MBPT2_R12 object
  R12IntEvalInfo(MBPT2_R12*);
  ~R12IntEvalInfo();

  void save_data_state(StateOut&);

  /** Sets whether to use dynamic load balancing in parallel MO transformations.
      Default is no */
  void set_dynamic(bool dynamic) { dynamic_ = dynamic; };
  /// Sets how frequently updates of progress are printed out. Default is 10%
  void set_print_percent(double print_percent) { print_percent_ = print_percent; };
  /// Set debug level. Default is 0.
  void set_debug_level(int debug) { debug_ = debug; };
  /** Sets the method of storing transformed MO integrals. Default depends on
      how the object was constructed. */
  void set_ints_method(const StoreMethod::type method) { ints_method_ = method; };
  /** Sets name of the file used to store transformed integrals.
      Default depends on how the object was constructed. */
  void set_ints_file(const std::string& filename) { ints_file_ = filename; };
  /** Sets the amount of memory to use for the calculation. Default is
      determined by DEFAULT_SC_MEMORY. */
  void set_memory(const size_t nbytes);
  /** Sets the ABS approach to be used (ABS or CABS).
      Default depends on how the object was constructed. */
  void set_absmethod(LinearR12::ABSMethod abs_method);

  Wavefunction* wfn() const { return wfn_; };
  Ref<Integral> integral() const { return refinfo()->ref()->integral(); };
  /// Returns the orbital basis set (OBS) object
  Ref<GaussianBasisSet> basis() const { return refinfo()->ref()->basis(); };
  /// Returns the virtuals basis set (VBS) obje19ct
  Ref<GaussianBasisSet> basis_vir() const { return bs_vir_; };
  /// Returns the resolution-of-the-identity basis set (RIBS) object
  Ref<GaussianBasisSet> basis_ri() const { return bs_ri_; };
  Ref<SCMatrixKit> matrixkit() const { return matrixkit_; };
  Ref<MemoryGrp> mem() const { return mem_;};
  Ref<MessageGrp> msg() const { return msg_;};
  Ref<ThreadGrp> thr() const { return thr_;};

  bool dynamic() const { return dynamic_; };
  double print_percent() const { return print_percent_; };
  int debug_level() const { return debug_; };
  const StoreMethod::type ints_method() const { return ints_method_; };
  const std::string& ints_file() const;
  const size_t memory() const { return memory_; };

  int nvir() const { return vir_->rank();};
  int nvir_act() const { return vir_act_->rank();};

  const Ref<LinearR12::CorrelationFactor>& corrfactor() const { return corrfactor_; }
  LinearR12::StandardApproximation stdapprox() const { return stdapprox_; }
  const Ref<LinearR12Ansatz>& ansatz() const { return ansatz_; }
  LinearR12::ABSMethod abs_method() const { return abs_method_; }
  bool gbc() const { return gbc_; }
  bool ebc() const { return ebc_; }
  bool spinadapted() const { return spinadapted_; }
  bool include_mp1() const { return include_mp1_; }
  unsigned int maxnabs() const { return maxnabs_; }
  bool ks_ebcfree() const { return ks_ebcfree_; }
  bool omit_P() const { return omit_P_; }
  bool safety_check() const { return safety_check_; }
  const LinearR12::PositiveDefiniteB& posdef_B() const { return posdef_B_; }

  /// Returns the MOIndexSpace object for all unoccupied MOs ordered by energy
  const Ref<MOIndexSpace>& vir() const { throw_if_spin_polarized(); return vir_; };
  /// Returns the MOIndexSpace object for all unoccupied MOs ordered by symmetry
  const Ref<MOIndexSpace>& vir_sb() const { throw_if_spin_polarized(); return vir_sb_; };
  /// Returns the MOIndexSpace object for the active unoccupied MOs
  const Ref<MOIndexSpace>& vir_act() const { throw_if_spin_polarized(); return vir_act_; };
  /// Returns the MOIndexSpace object for all unoccupied MOs ordered by energy
  const Ref<MOIndexSpace>& vir(const SpinCase1& S) const { return vir_spaces_[S].vir_; };
  /// Returns the MOIndexSpace object for all unoccupied MOs ordered by symmetry
  const Ref<MOIndexSpace>& vir_sb(const SpinCase1& S) const { return vir_spaces_[S].vir_sb_; };
  /// Returns the MOIndexSpace object for the active unoccupied MOs
  const Ref<MOIndexSpace>& vir_act(const SpinCase1& S) const { return vir_spaces_[S].vir_act_; };

  /// Cheating! fock() is not available yet standalone, thus these spaces must be modified after canonicalization
  void vir(const SpinCase1& S, const Ref<MOIndexSpace>& space);
  void vir_sb(const SpinCase1& S, const Ref<MOIndexSpace>& space);
  void vir_act(const SpinCase1& S, const Ref<MOIndexSpace>& space);
  
  /// Returns the MOIndexSpace object for ABS
  const Ref<MOIndexSpace>& abs_space() const { return abs_space_; };
  /// Returns the MOIndexSpace object for RI-BS: approximates the identity
  const Ref<MOIndexSpace>& ribs_space() const { return ribs_space_; };
  /** Returns the MOIndexSpace object for RI-BS:
      if CABS/CABS+ -- approximates the complement to OBS,
      if ABS/ABS+   -- null
  */
  const Ref<MOIndexSpace>& ribs_space(const SpinCase1& S) const;
  /// Returns the MOIntsTransformFactory object
  const Ref<MOIntsTransformFactory>& tfactory() const { return tfactory_; };
  /// Return the SingleRefInfo object
  const Ref<SingleRefInfo>& refinfo() const;
  
  /** Compute span of bs and create corresponding mospace referred to by name. Number
      linear dependencies is returned in nlindep */
  static Ref<MOIndexSpace> orthogonalize(const std::string& id, const std::string& name, const Ref<GaussianBasisSet>& bs,
                                         const Ref<Integral>& integral, OverlapOrthog::OrthogMethod orthog_method, double lindep_tol,
                                         int& nlindep);

  /** Project space1 on space2. This routine computes X2 such that C1.S12.X2 = I,
      where I is identity matrix, C1 is space1, and X2 spans
      subspace of space2. X2 is returned. */
  static Ref<MOIndexSpace> gen_project(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                       const std::string& id, const std::string& name, double lindep_tol);
  /** Compute subspace X2 of space2 which is orthogonal complement to space1, i.e.,
      C1.S12.X2=0, where 0 is the null matrix.
  */
  static Ref<MOIndexSpace> orthog_comp(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                const std::string& id, const std::string& name, double lindep_tol);
                                       
  /// Compute overlap matrices in the basis of space1 and space2
  static void compute_overlap_ints(const Ref<MOIndexSpace>& space1,
                            const Ref<MOIndexSpace>& space2,
                            RefSCMatrix& S);
  /// Compute electric dipole and quadrupole moment matrices in the basis of space1 and space2
  static void compute_multipole_ints(const Ref<MOIndexSpace>& space1,
                              const Ref<MOIndexSpace>& space2,
                              RefSCMatrix& MX,
                              RefSCMatrix& MY,
                              RefSCMatrix& MZ,
                              RefSCMatrix& MXX,
                              RefSCMatrix& MYY,
                              RefSCMatrix& MZZ);
			      
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


