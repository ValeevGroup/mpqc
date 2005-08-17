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

#ifndef _chemistry_qc_mbptr12_mp2r12energy_h
#define _chemistry_qc_mbptr12_mp2r12energy_h

#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/linearr12.h>
//#include <chemistry/qc/mbptr12/vxb_eval.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/twobodygrid.h>

namespace sc {

  /** Class MP2R12Energy is the object that computes and maintains MP2-R12 energies */

class MP2R12Energy : virtual public SavableState {
  private:
  Ref<R12IntEval> r12eval_;
  LinearR12::StandardApproximation stdapprox_;
  bool ebc_;
  int debug_;
  bool evaluated_;
  
  RefSCVector er12_aa_, er12_ab_, emp2r12_aa_, emp2r12_ab_;
  RefSCVector ef12_[NSpinCases2], emp2f12_[NSpinCases2];
  // The coefficients are stored ij by kl, where kl is the r12-multiplied pair
  RefSCMatrix Caa_, Cab_;
  RefSCMatrix C_[NSpinCases2];

  double emp2tot_aa_() const;
  double emp2tot_ab_() const;
  double er12tot_aa_();
  double er12tot_ab_();
  double emp2tot_(SpinCase2 S) const;
  double ef12tot_(SpinCase2 S) const;

  // Initialize SCVectors and SCMatrices
  void init_();
  // new compute function to replace old compute() in the future
  void compute_new_();

  // Computes values of all 2-body products from
  // space1 and space2 if electron 1 is at r1 and
  // electron 2 is at r2. equiv specifies whether electrons
  // are equivalent (same spin) or not
  RefSCVector compute_2body_values_(bool equiv, const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                    const SCVector3& r1, const SCVector3& r2) const;

public:

  MP2R12Energy(StateIn&);
  MP2R12Energy(Ref<R12IntEval>& r12eval, LinearR12::StandardApproximation stdapp, int debug);
  ~MP2R12Energy();

  void save_data_state(StateOut&);
  void obsolete();
  void print(std::ostream&o=ExEnv::out0()) const;
  void print_pair_energies(bool spinadapted, std::ostream&so=ExEnv::out0());

  Ref<R12IntEval> r12eval() const;
  LinearR12::StandardApproximation stdapp() const;
  /** Returns whether Generalized Brillouin Condition (GBC) was used in evaluation of
      the MP2-R12 intermediates */
  const bool gbc() const;
  /** Returns whether Extended Brillouin Condition (EBC) was used in evaluation of
      the MP2-R12 intermediates and the MP2-R12 energy */
  const bool ebc() const;
  void set_debug(int debug);
  int get_debug() const;
  
  /// Computes the first-order R12 wave function and MP2-R12 energy
  void compute();
  /** Computes the value of the alpha-alpha pair function ij
      when electrons 1 and 2 reside at r1 and r2 */
  double compute_pair_function_aa(int ij, const SCVector3& r1, const SCVector3& r2);
  /** Computes the value of the alpha-beta pair function ij
      when electrons 1 and 2 reside at r1 and r2 */
  double compute_pair_function_ab(int ij, const SCVector3& r1, const SCVector3& r2);
  /** Computes values of the alpha-alpha pair function ij on tbgrid */
  void compute_pair_function_aa(int ij, const Ref<TwoBodyGrid>& tbgrid);
  /** Computes values of the alpha-beta pair function ij on tbgrid */
  void compute_pair_function_ab(int ij, const Ref<TwoBodyGrid>& tbgrid);

  /// Returns the vector of MP2 alpha-alpha pair energies
  RefSCVector emp2_aa() const;
  /// Returns the vector of MP2 alpha-beta pair energies
  RefSCVector emp2_ab() const;
  /// Returns the vector of R12 corrections to MP2-R12 alpha-alpha pair energies
  RefSCVector er12_aa() const;
  /// Returns the vector of R12 correction to MP2-R12 alpha-beta pair energies
  RefSCVector er12_ab() const;
  /// Returns the vector of MP2-R12 alpha-alpha pair energies
  RefSCVector emp2r12_aa() const;
  /// Returns the vector of MP2-R12 alpha-beta pair energies
  RefSCVector emp2r12_ab() const;
  /// Returns the vector of second-order pair energies of spin case S
  RefSCVector emp2(SpinCase2 S) const;
  /// Returns the vector of F12 corrections to second-order pair energies of spin case S
  RefSCVector ef12(SpinCase2 S) const;
  /// Returns total MP2-F12 correlation energy
  double energy();

  /** Returns the matrix of amplitudes of
      alpha-alpha r12-multiplied occupied orbital pairs in the first-order
      pair function
  */
  RefSCMatrix C_aa();
  /** Returns the matrix of amplitudes of
      alpha-beta r12-multiplied occupied orbital pairs in the first-order
      pair function
  */
  RefSCMatrix C_ab();
  /** Returns the matrix of amplitudes of
      alpha-alpha virtuals orbital pairs in the first-order
      pair function
  */
  RefSCMatrix T2_aa();
  /** Returns the matrix of amplitudes of
      alpha-beta virtuals orbital pairs in the first-order
      pair function
  */
  RefSCMatrix T2_ab();

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


