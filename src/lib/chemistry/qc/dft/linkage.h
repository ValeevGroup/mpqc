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
#include <chemistry/qc/dft/solvent.h>

namespace sc {

ForceLink<RadialAngularIntegrator> dft_force_link_a_;
ForceLink<NElFunctional> dft_force_link_b_;
ForceLink<XalphaFunctional> dft_force_link_c_;
ForceLink<SlaterXFunctional> dft_force_link_d_;
ForceLink<Becke88XFunctional> dft_force_link_e_;
ForceLink<LYPCFunctional> dft_force_link_f_;
ForceLink<CLKS> dft_force_link_h_;
ForceLink<UKS> dft_force_link_i_;
ForceLink<VWN5LCFunctional> dft_force_link_j_;
ForceLink<VWN3LCFunctional> dft_force_link_k_;
ForceLink<PW92LCFunctional> dft_force_link_l_;
ForceLink<PBEXFunctional> dft_force_link_m_;
ForceLink<PBECFunctional> dft_force_link_n_;
ForceLink<P86CFunctional> dft_force_link_o_;
ForceLink<PW91XFunctional> dft_force_link_p_;
ForceLink<PW86XFunctional> dft_force_link_q_;
ForceLink<PZ81LCFunctional> dft_force_link_r_;
ForceLink<G96XFunctional> dft_force_link_s_;
ForceLink<VWN1LCFunctional> dft_force_link_t_;
ForceLink<VWN2LCFunctional> dft_force_link_u_;
ForceLink<VWN4LCFunctional> dft_force_link_v_;
ForceLink<PW91CFunctional> dft_force_link_w_;
ForceLink<HSOSKS> dft_force_link_x_;
ForceLink<VWNLCFunctional> dft_force_link_y_;
ForceLink<NewP86CFunctional> dft_force_link_z_;
ForceLink<AM05Functional> dft_force_link_aa_;
ForceLink<BEMSolventH> dft_force_link_ab_;

}

#endif
