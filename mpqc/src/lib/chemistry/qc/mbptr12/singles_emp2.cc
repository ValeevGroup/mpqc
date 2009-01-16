//
// singles_emp2.cc
//
// Copyright (C) 2007 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_singles_emp2_()
{
  if (evaluated_)
    return;
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();

  Timer tim("singles MP2 energy");
  ExEnv::out0() << endl << indent
	       << "Entered singles MP2 energy evaluator" << endl;
  ExEnv::out0() << incindent;

  int me = msg->me();
  int nproc = msg->n();
  
  emp2_singles_ = 0.0;
  for(int s=0; s<nspincases1(); s++) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    
    Ref<OrbitalSpace> occ_act = r12info_->refinfo()->occ_act(spin);
    Ref<OrbitalSpace> vir_act = r12info_->vir_act(spin);
    RefSCMatrix Fia = fock(occ_act,vir_act,spin);

    const int ni = occ_act->rank();
    const int na = vir_act->rank();
    const RefDiagSCMatrix& ievals = occ_act->evals();
    const RefDiagSCMatrix& aevals = vir_act->evals();

    for(int i=0; i<ni; i++) {
      for(int a=0; a<na; a++) {
        const double fia = Fia.get_element(i,a);
        emp2_singles_ -= fia*fia/(-ievals(i)+aevals(a));
      }
    }
  }
  
  ExEnv::out0() << decindent;
  ExEnv::out0() << endl << "Exited singles MP2 energy evaluator" << endl;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
