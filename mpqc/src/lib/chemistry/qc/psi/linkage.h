//
// linkage.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_psi_linkage_h
#define _chemistry_qc_psi_linkage_h

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psicc.h>
#include <chemistry/qc/psi/psici.h>
#include <chemistry/qc/psi/psicc_pt2r12.h>
#include <chemistry/qc/psi/psici_pt2r12.h>
#include <math/optimize/qnewton.h>

namespace sc {

static ForceLink<PsiCLHF> psi_force_link_a_;
static ForceLink<PsiHSOSHF> psi_force_link_b_;
static ForceLink<PsiUHF> psi_force_link_c_;
static ForceLink<PsiCCSD> psi_force_link_d_;
static ForceLink<PsiCCSD_T> psi_force_link_e_;
static ForceLink<PsiCCSD_PT2R12> psi_force_link_f_;
static ForceLink<PsiCCSD_PT2R12T> psi_force_link_g_;
static ForceLink<PsiCI> psi_force_link_h_;
static ForceLink<PsiCI_PT2R12> psi_force_link_i_;

}

#endif
