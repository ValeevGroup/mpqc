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
#include <chemistry/qc/mbptr12/vxb_eval.h>

namespace sc {

  /** Class MP2R12Energy is the object that computes and maintains MP2-R12 energies */

class MP2R12Energy : virtual public SavableState {

  Ref<R12IntEval> r12eval_;
  LinearR12::StandardApproximation stdapprox_;
  int debug_;
  bool evaluated_;
  
  RefSCVector er12_aa_, er12_ab_, emp2r12_aa_, emp2r12_ab_;

  double emp2tot_aa_() const;
  double emp2tot_ab_() const;
  double er12tot_aa_();
  double er12tot_ab_();

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
  void set_debug(int debug);
  int get_debug() const;
  
  RefSCDimension dim_aa() const;
  RefSCDimension dim_ab() const;
  RefSCDimension dim_s() const;
  RefSCDimension dim_t() const;

  void compute();
  
  RefSCVector emp2_aa() const;
  RefSCVector emp2_ab() const;
  RefSCVector er12_aa() const;
  RefSCVector er12_ab() const;
  RefSCVector emp2r12_aa() const;
  RefSCVector emp2r12_ab() const;

  /// Total MP2-R12 correlation energy
  double energy();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


