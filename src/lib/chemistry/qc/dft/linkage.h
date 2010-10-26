//
// linkage.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#ifndef _chemistry_qc_dft_linkage_h
#define _chemistry_qc_dft_linkage_h

#include <chemistry/qc/dft/clks.h>
#include <chemistry/qc/dft/uks.h>
#include <chemistry/qc/dft/hsosks.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/dft/am05.h>

namespace sc {

static ForceLink<RadialAngularIntegrator> dft_force_link_a_;
static ForceLink<NElFunctional> dft_force_link_b_;
static ForceLink<XalphaFunctional> dft_force_link_c_;
static ForceLink<SlaterXFunctional> dft_force_link_d_;
static ForceLink<Becke88XFunctional> dft_force_link_e_;
static ForceLink<LYPCFunctional> dft_force_link_f_;
static ForceLink<CLKS> dft_force_link_h_;
static ForceLink<UKS> dft_force_link_i_;
static ForceLink<VWN5LCFunctional> dft_force_link_j_;
static ForceLink<VWN3LCFunctional> dft_force_link_k_;
static ForceLink<PW92LCFunctional> dft_force_link_l_;
static ForceLink<PBEXFunctional> dft_force_link_m_;
static ForceLink<PBECFunctional> dft_force_link_n_;
static ForceLink<P86CFunctional> dft_force_link_o_;
static ForceLink<PW91XFunctional> dft_force_link_p_;
static ForceLink<PW86XFunctional> dft_force_link_q_;
static ForceLink<PZ81LCFunctional> dft_force_link_r_;
static ForceLink<G96XFunctional> dft_force_link_s_;
static ForceLink<VWN1LCFunctional> dft_force_link_t_;
static ForceLink<VWN2LCFunctional> dft_force_link_u_;
static ForceLink<VWN4LCFunctional> dft_force_link_v_;
static ForceLink<PW91CFunctional> dft_force_link_w_;
static ForceLink<HSOSKS> dft_force_link_x_;
static ForceLink<VWNLCFunctional> dft_force_link_y_;
static ForceLink<NewP86CFunctional> dft_force_link_z_;
static ForceLink<AM05Functional> dft_force_link_aa_;

}

#endif
