//
// linkage.h
//
// Copyright (C) 2001 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_linkage_h
#define _chemistry_qc_mbptr12_linkage_h

#include <util/class/class.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/pt2r12.h>
#include <chemistry/qc/mbptr12/molcas_pt2r12.h>

#include <chemistry/qc/mbpt/linkage.h>

namespace sc {

ForceLink<MBPT2_R12> mbptr12_force_link_a_;
ForceLink<PT2R12> mbptr12_force_link_b_;
ForceLink<SpinOrbitalPT2R12> mbptr12_force_link_c_;
ForceLink<MolcasPT2R12> mbptr12_force_link_d_;

}

#endif
