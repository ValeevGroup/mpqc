//
// compute_amps.cc
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
#include <string.h>
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
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_T2_(RefSCMatrix& T2,
                        const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2,
                        const Ref<MOIndexSpace>& space3,
                        const Ref<MOIndexSpace>& space4,
                        bool antisymmetrize,
                        const Ref<TwoBodyMOIntsTransform>& transform)
{
  tim_enter("T2 amplitudes");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|T2|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
	       << "Entered MP2 T2 amplitude (" << label << ") evaluator" << endl;
  ExEnv::out0() << incindent;
  
  typedef std::vector< Ref<TwoBodyMOIntsTransform> > tformvec;
  tformvec tform;
  if (transform.nonnull())
    tform.push_back(transform);
  compute_tbint_tensor<ManyBodyTensors::ERI_to_T2,false,false>(
    T2, corrfactor()->tbint_type_eri(),
    space1, space2,
    space3, space4,
    antisymmetrize,
    tform,
    std::vector< Ref<TwoBodyIntDescr> >(1,new TwoBodyIntDescrERI(r12info()->integral()))
  );
  
  ExEnv::out0() << decindent << indent << "Exited MP2 T2 amplitude (" << label << ") evaluator" << endl;
  tim_exit("T2 amplitudes");
}


void
R12IntEval::compute_F12_(RefSCMatrix& F12,
                        const Ref<MOIndexSpace>& space1,
                        const Ref<MOIndexSpace>& space2,
                        const Ref<MOIndexSpace>& space3,
                        const Ref<MOIndexSpace>& space4,
                        const std::vector< Ref<TwoBodyMOIntsTransform> >& transvec,
                        const std::vector< Ref<TwoBodyIntDescr> >& descrvec)
{
  tim_enter("F12 amplitudes");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|F12|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
                << "Entered F12 amplitude (" << label << ") evaluator" << endl
                << incindent;
  
  compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
    F12, corrfactor()->tbint_type_f12(),
    space1, space2,
    space3, space4,
    false,
    transvec, descrvec
  );

  ExEnv::out0() << decindent << indent << "Exited F12 amplitude (" << label << ") evaluator" << endl;
  tim_exit("F12 amplitudes");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
