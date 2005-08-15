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

#define USE_SINGLEREFINFO 1

#include <util/misc/string.h>
#include <util/ref/ref.h>
#include <math/scmat/abstract.h>
#include <util/group/memory.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/moindexspace.h>
#include <chemistry/qc/mbptr12/transform_factory.h>
#if USE_SINGLEREFINFO
#  include <chemistry/qc/mbptr12/singlerefinfo.h>
#endif

#ifndef _chemistry_qc_mbptr12_vxbevalinfo_h
#define _chemistry_qc_mbptr12_vxbevalinfo_h

namespace sc {

class MBPT2_R12;

  /** Class R12IntEvalInfo contains information necessary for R12 intermediate
      evaluators */

class R12IntEvalInfo : virtual public SavableState {

public:

  /// Describes the method of storing transformed MO integrals. See MBPT2_R12.
  enum StoreMethod { mem_posix = 0, posix = 1, mem_mpi = 2, mpi = 3, mem_only = 4 };

private:

  Wavefunction* wfn_;     // Wave function that owns this
#if !USE_SINGLEREFINFO
  Ref<SCF> ref_;
  Ref<Integral> integral_;
  Ref<GaussianBasisSet> bs_;
#endif
  Ref<GaussianBasisSet> bs_aux_;
  Ref<GaussianBasisSet> bs_vir_;
  Ref<GaussianBasisSet> bs_ri_;
  Ref<SCMatrixKit> matrixkit_;
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;

#if !USE_SINGLEREFINFO
  int ndocc_;
  int nsocc_;
  int nfzc_;
  int nfzv_;
#endif

  size_t memory_;
  bool dynamic_;
  double print_percent_;
  int debug_;
  StoreMethod ints_method_;
  char* ints_file_;
  LinearR12::ABSMethod abs_method_;

  int nlindep_aux_;
  int nlindep_vir_;
  int nlindep_ri_;
  
  /// This space depends only on the orbital basis set
  Ref<MOIndexSpace> abs_space_;  // ABS space
  Ref<MOIndexSpace> ribs_space_; // RIBS basis  

#if !USE_SINGLEREFINFO
  Ref<MOIndexSpace> mo_space_;   // symblocked MO space
  Ref<MOIndexSpace> obs_space_;  // energy-sorted MO space
  Ref<MOIndexSpace> act_occ_space_;
  Ref<MOIndexSpace> occ_space_;
  Ref<MOIndexSpace> occ_space_symblk_;
#endif
  Ref<MOIndexSpace> act_vir_space_;
  Ref<MOIndexSpace> vir_space_;
  Ref<MOIndexSpace> vir_space_symblk_;
  
#if USE_SINGLEREFINFO
  /// Initializes all spaces that relate to the reference determinant
  Ref<SingleRefInfo> refinfo_;
#endif
  
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
#if !USE_SINGLEREFINFO
  // Construct eigenvector and eigenvalues sorted by energy
  void eigen_();
#endif
  // Construct orthog_aux_
  void construct_orthog_aux_();
  // Construct orthog_vir_
  void construct_orthog_vir_();
  // Construct orthog_ri_
  void construct_orthog_ri_();

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
  void set_ints_method(const StoreMethod method) { ints_method_ = method; };
  /** Sets name of the file used to store transformed integrals.
      Default depends on how the object was constructed. */
  void set_ints_file(const char* filename) { ints_file_ = strdup(filename); };
  /** Sets the amount of memory to use for the calculation. Default is
      determined by DEFAULT_SC_MEMORY. */
  void set_memory(const size_t nbytes);
  /** Sets the ABS approach to be used (ABS or CABS).
      Default depends on how the object was constructed. */
  void set_absmethod(LinearR12::ABSMethod abs_method);

  Wavefunction* wfn() const { return wfn_; };
#if !USE_SINGLEREFINFO
  Ref<SCF> ref() const { return ref_; };
  Ref<Integral> integral() const { return integral_; };
  /// Returns the orbital basis set (OBS) object
  Ref<GaussianBasisSet> basis() const { return bs_; };
#else
  Ref<Integral> integral() const { return refinfo()->ref()->integral(); };
  /// Returns the orbital basis set (OBS) object
  Ref<GaussianBasisSet> basis() const { return refinfo()->ref()->basis(); };
#endif
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
  const StoreMethod ints_method() const { return ints_method_; };
  char* ints_file() const;
  const size_t memory() const { return memory_; };

#if !USE_SINGLEREFINFO
  int ndocc() const { return ndocc_;};
  int nsocc() const { return nsocc_;};
  int nfzc() const { return nfzc_;};
  int nfzv() const { return nfzv_;};
  int ndocc_act() const { return ndocc() - nfzc();};
  int nalpha() const { return ndocc()+ nsocc(); }
  int nbeta() const { return ndocc(); }
  int nalpha_occ() const { return ndocc_act() + nsocc(); }
  int nbeta_occ() const { return ndocc_act(); }
#endif
  int nvir() const { return vir_space_->rank();};
  int nvir_act() const { return act_vir_space_->rank();};

  LinearR12::ABSMethod abs_method() const { return abs_method_; };

#if !USE_SINGLEREFINFO
  /// Returns the MOIndexSpace object for symmetry-blocked MOs in OBS
  const Ref<MOIndexSpace>& mo_space() const { return mo_space_; };
  /// Returns the MOIndexSpace object for energy-sorted MOs in OBS
  const Ref<MOIndexSpace>& obs_space() const { return obs_space_; };
  /// Returns the MOIndexSpace object for the active occupied MOs
  const Ref<MOIndexSpace>& act_occ_space() const { return act_occ_space_; };
  /// Returns the MOIndexSpace object for all occupied MOs sorted by energy
  const Ref<MOIndexSpace>& occ_space() const { return occ_space_; };
  /// Returns the MOIndexSpace object for all occupied MOs symmetry-blocked
  const Ref<MOIndexSpace>& occ_space_symblk() const { return occ_space_symblk_; };
#endif
  /// Returns the MOIndexSpace object for all unoccupied MOs ordered by energy
  const Ref<MOIndexSpace>& vir_space() const { return vir_space_; };
  /// Returns the MOIndexSpace object for all unoccupied MOs ordered by symmetry
  const Ref<MOIndexSpace>& vir_space_symblk() const { return vir_space_symblk_; };
  /// Returns the MOIndexSpace object for the active unoccupied MOs
  const Ref<MOIndexSpace>& act_vir_space() const { return act_vir_space_; };
  /// Returns the MOIndexSpace object for ABS
  const Ref<MOIndexSpace>& abs_space() const { return abs_space_; };
  /// Returns the MOIndexSpace object for RI-BS
  const Ref<MOIndexSpace>& ribs_space() const { return ribs_space_; };
  /// Returns the MOIntsTransformFactory object
  const Ref<MOIntsTransformFactory>& tfactory() const { return tfactory_; };
#if USE_SINGLEREFINFO
  /// Return the SingleRefInfo object
  const Ref<SingleRefInfo>& refinfo() const;
#endif
  
  /** Compute span of bs and create corresponding mospace referred to by name. Number
      linear dependencies is returned in nlindep */
  static Ref<MOIndexSpace> orthogonalize(const std::string& id, const std::string& name, const Ref<GaussianBasisSet>& bs,
                                  OverlapOrthog::OrthogMethod orthog_method, double lindep_tol,
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


