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

#ifndef _chemistry_qc_mbptr12_vxbevalinfo_h
#define _chemistry_qc_mbptr12_vxbevalinfo_h

#include <util/ref/ref.h>
#include <math/scmat/abstract.h>
#include <util/group/memory.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/scf/scf.h>

namespace sc {

class MBPT2_R12;

  /** Class R12IntEvalInfo contains information necessary for R12 intermediate
      evaluators */

class R12IntEvalInfo : virtual public SavableState {

public:

  enum StoreMethod { mem_posix = 0, posix = 1, mem_mpi = 2, mpi = 3, mem_only = 4 };

private:

  MolecularEnergy* mole_;     // MolecularEnergy that owns this
  Ref<SCF> ref_;
  Ref<Integral> integral_;
  Ref<GaussianBasisSet> bs_;
  Ref<GaussianBasisSet> bs_aux_;
  Ref<SCMatrixKit> matrixkit_;
  Ref<MessageGrp> msg_;
  Ref<MemoryGrp> mem_;
  Ref<ThreadGrp> thr_;

  int nocc_;
  int nocc_act_;
  int nfzc_;
  int nfzv_;
  int noso_;

  size_t memory_;
  bool dynamic_;
  int debug_;
  StoreMethod ints_method_;
  char* ints_file_;

  RefSCMatrix scf_vec_;
  RefDiagSCMatrix evals_;
  RefDiagSCMatrix occs_;
  int* orbsym_;

  void eigen_(RefDiagSCMatrix& evals, RefSCMatrix& scf_vec, RefDiagSCMatrix& occs, int*& orbsym);

public:
  R12IntEvalInfo(StateIn&);
  R12IntEvalInfo(MBPT2_R12*);
  ~R12IntEvalInfo();

  void save_data_state(StateOut&);

  void set_dynamic(bool dynamic) { dynamic_ = dynamic; };
  void set_debug_level(int debug) { debug_ = debug; };
  void set_ints_method(const StoreMethod method) { ints_method_ = method; };
  void set_ints_file(const char* filename) { ints_file_ = strdup(filename); };
  void set_memory(size_t nbytes) { if (nbytes >= 0) memory_ = nbytes; };

  MolecularEnergy* mole() const { return mole_; };
  Ref<SCF> ref() const { return ref_; };
  Ref<Integral> integral() const { return integral_; };
  Ref<GaussianBasisSet> basis() const { return bs_; };
  Ref<GaussianBasisSet> basis_aux() const { return bs_aux_; };
  Ref<SCMatrixKit> matrixkit() const { return matrixkit_; };
  Ref<MemoryGrp> mem() const { return mem_;};
  Ref<MessageGrp> msg() const { return msg_;};
  Ref<ThreadGrp> thr() const { return thr_;};

  bool dynamic() const { return dynamic_; };
  int debug_level() const { return debug_; };
  const StoreMethod ints_method() const { return ints_method_; };
  const char* ints_file() const;
  const size_t memory() const { return memory_; };

  const int nocc() const { return nocc_;};
  const int nocc_act() const { return nocc_act_;};
  const int noso() const { return noso_;};
  const int nfzc() const { return nfzc_;};
  const int nfzv() const { return nfzv_;};

  RefSCMatrix scf_vec() const { return scf_vec_; };
  RefDiagSCMatrix evals() const { return evals_; };
  int *orbsym() const { return orbsym_; };
  
  /// Compute dipole and quadrupole moment matrices in active MO basis
  void compute_multipole_ints(RefSymmSCMatrix& MX,
			      RefSymmSCMatrix& MY,
			      RefSymmSCMatrix& MZ,
			      RefSymmSCMatrix& MXX,
			      RefSymmSCMatrix& MYY,
			      RefSymmSCMatrix& MZZ);
			      
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


