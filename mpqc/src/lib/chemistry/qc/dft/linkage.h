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

#ifndef __PIC__

#include <chemistry/qc/dft/clks.h>
#include <chemistry/qc/dft/uks.h>
#include <chemistry/qc/dft/hsosks.h>
#include <chemistry/qc/dft/hfacm.h>
#include <chemistry/qc/dft/integrator.h>
#include <chemistry/qc/dft/functional.h>

const ClassDesc &dft_force_link_a_ = Murray93Integrator::class_desc_;
const ClassDesc &dft_force_link_b_ = NElFunctional::class_desc_;
const ClassDesc &dft_force_link_c_ = XalphaFunctional::class_desc_;
const ClassDesc &dft_force_link_d_ = SlaterXFunctional::class_desc_;
const ClassDesc &dft_force_link_e_ = Becke88XFunctional::class_desc_;
const ClassDesc &dft_force_link_f_ = LYPCFunctional::class_desc_;
const ClassDesc &dft_force_link_g_ = HFACM::class_desc_;
const ClassDesc &dft_force_link_h_ = CLKS::class_desc_;
const ClassDesc &dft_force_link_i_ = UKS::class_desc_;
const ClassDesc &dft_force_link_j_ = VWN5LCFunctional::class_desc_;
const ClassDesc &dft_force_link_k_ = VWN3LCFunctional::class_desc_;
const ClassDesc &dft_force_link_l_ = PW92LCFunctional::class_desc_;
const ClassDesc &dft_force_link_m_ = PBEXFunctional::class_desc_;
const ClassDesc &dft_force_link_n_ = PBECFunctional::class_desc_;
const ClassDesc &dft_force_link_o_ = P86CFunctional::class_desc_;
const ClassDesc &dft_force_link_p_ = PW91XFunctional::class_desc_;
const ClassDesc &dft_force_link_q_ = PW86XFunctional::class_desc_;
const ClassDesc &dft_force_link_r_ = PZ81LCFunctional::class_desc_;
const ClassDesc &dft_force_link_s_ = G96XFunctional::class_desc_;
const ClassDesc &dft_force_link_t_ = VWN1LCFunctional::class_desc_;
const ClassDesc &dft_force_link_u_ = VWN2LCFunctional::class_desc_;
const ClassDesc &dft_force_link_v_ = VWN4LCFunctional::class_desc_;
const ClassDesc &dft_force_link_w_ = PW91CFunctional::class_desc_;
const ClassDesc &dft_force_link_x_ = HSOSKS::class_desc_;
// const ClassDesc &dft_force_link_y_ = EMLIntegrator::class_desc_;

#endif /* __PIC__ */

#endif
