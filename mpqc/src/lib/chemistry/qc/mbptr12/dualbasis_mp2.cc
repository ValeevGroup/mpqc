//
// dualbasis_mp2.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
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
R12IntEval::compute_dualEmp1_()
{
  if (evaluated_)
    return;
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();

  tim_enter("dual-basis MP1 energy");
  ExEnv::out0() << endl << indent
	       << "Entered dual-basis MP1 energy evaluator" << endl;
  ExEnv::out0() << incindent;

  int me = msg->me();
  int nproc = msg->n();
  
  // Compute act.occ./aux.virt. Fock matrix
  form_canonvir_space_();
#if USE_SINGLEREFINFO
  Ref<MOIndexSpace> occ_space = r12info_->refinfo()->docc();
  const double eref = r12info_->refinfo()->ref()->energy();
#else
  Ref<MOIndexSpace> occ_space = r12info_->occ_space();
  const double eref = r12info_->ref()->energy();
#endif
  RefSCMatrix F_aocc_canonvir = fock_(occ_space,occ_space,canonvir_space_);

  int nocc = occ_space->rank();
  int ncanonvir = canonvir_space_->rank();
  RefDiagSCMatrix occ_evals = occ_space->evals();
  RefDiagSCMatrix canonvir_evals = canonvir_space_->evals();

  double emp1 = 0.0;
  for(int i=0; i<nocc; i++) {
    for(int a=0; a<ncanonvir; a++) {
      const double Fia = F_aocc_canonvir.get_element(i,a);
      emp1 += Fia*Fia/(-occ_evals(i)+canonvir_evals(a));
    }
  }
  ExEnv::out0() << indent << "MP1 energy correction to HF energy [au] :   "
                << 2.0*emp1 << endl;
  ExEnv::out0() << indent << "HF energy estimated in new basis [au]   :   "
                << eref - 2.0*emp1 << endl;

  ExEnv::out0() << decindent;
  ExEnv::out0() << endl << "Exited dual-basis MP1 energy evaluator" << endl;

  tim_exit("dual-basis MP1 energy");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
