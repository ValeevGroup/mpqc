//
// compute_a_gebc_vbs.cc
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
R12IntEval::contrib_to_VXB_gebc_vbsneqobs_()
{
  if (evaluated_)
    return;

  /*form_canonvir_space_();
  contrib_to_VXB_a_symm_("(im|jn)",r12info_->occ_space());
  contrib_to_VXB_a_symm_("(ia|jb)",canonvir_space);
  contrib_to_VXB_a_asymm_("(im|ja)",r12info_->occ_space(),canonvir_space_);
  if (r12info_->basis_vir() != r12info_->basis_ri())
    contrib_to_VXB_a_asymm_("(im|jy)",r12info_->occ_space(),r12info_->ribs_space());
  */
  contrib_to_VXB_a_symm_("(im|jn)");
  contrib_to_VXB_a_symm_("(ia|jb)");
  contrib_to_VXB_a_asymm_("(im|ja)");
  if (r12info_->basis_vir() != r12info_->basis_ri())
    contrib_to_VXB_a_asymm_("(im|jy)");

  
  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
