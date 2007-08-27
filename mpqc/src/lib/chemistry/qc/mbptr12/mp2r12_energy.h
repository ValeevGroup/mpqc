//
// mp2r12_energy.h
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

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/linearr12.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>

#ifndef _chemistry_qc_mbptr12_mp2r12energy_h
#define _chemistry_qc_mbptr12_mp2r12energy_h

#define MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION 1

namespace sc {

  /** Class MP2R12Energy is the object that computes and maintains MP2-R12 energies */
class MP2R12Energy : virtual public SavableState {
  private:
  Ref<R12IntEval> r12eval_;
  LinearR12::StandardApproximation stdapprox_;
  bool ebc_;
  int debug_;
  bool evaluated_;
  
  RefSCVector ef12_[NSpinCases2], emp2f12_[NSpinCases2];
  // The coefficients are stored xy by ij, where xy is the geminal-multiplied pair
  RefSCMatrix C_[NSpinCases2];

  double emp2f12tot(SpinCase2 S) const;
  double ef12tot(SpinCase2 S) const;

  // Initialize SCVectors and SCMatrices
  void init_();
  /** Computes values of all 2-body products from
      space1 and space2 if electron 1 is at r1 and
      electron 2 is at r2. equiv specifies whether electrons
      are equivalent (same spin) or not */
  RefSCVector compute_2body_values_(bool equiv, const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                    const SCVector3& r1, const SCVector3& r2) const;

public:

  MP2R12Energy(StateIn&);
  MP2R12Energy(Ref<R12IntEval>& r12eval, LinearR12::StandardApproximation stdapp, int debug);
  ~MP2R12Energy();

  void save_data_state(StateOut&);
  void obsolete();
  void print(std::ostream&o=ExEnv::out0()) const;

  Ref<R12IntEval> r12eval() const;
  LinearR12::StandardApproximation stdapp() const;
  /** Returns whether Generalized Brillouin Condition (GBC) was used in evaluation of
      the MP2-R12 intermediates */
  bool gbc() const;
  /** Returns whether Extended Brillouin Condition (EBC) was used in evaluation of
      the MP2-R12 intermediates and the MP2-R12 energy */
  bool ebc() const;
  void set_debug(int debug);
  int get_debug() const;
  
  /// Computes the first-order R12 wave function and MP2-R12 energy
  void compute();
  // Print pair energies nicely to so
  void print_pair_energies(bool spinadapted, std::ostream&so=ExEnv::out0());

#if MP2R12ENERGY_CAN_COMPUTE_PAIRFUNCTION
  /** Computes values of pair function ij on tbgrid */
  void compute_pair_function(unsigned int i, unsigned int j, SpinCase2 spincase2,
                             const Ref<TwoBodyGrid>& tbgrid);
#endif

  /// Returns the vector of second-order pair energies of spin case S
  const RefSCVector& emp2f12(SpinCase2 S) const;
  /// Returns the vector of F12 corrections to second-order pair energies of spin case S
  const RefSCVector& ef12(SpinCase2 S) const;
  /// Returns total MP2-F12 correlation energy
  double energy();

  /** Returns the matrix of first-order amplitudes of r12-multiplied occupied orbital pairs.
  */
  RefSCMatrix C(SpinCase2 S);
  /** Returns the matrix of first-order amplitudes of conventional orbital pairs.
  */
  RefSCMatrix T2(SpinCase2 S);

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


