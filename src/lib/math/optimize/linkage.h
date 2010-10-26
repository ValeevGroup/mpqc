//
// linkage.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _math_optimize_linkage_h
#define _math_optimize_linkage_h

#include <math/optimize/qnewton.h>
#include <math/optimize/newton.h>
#include <math/optimize/gdiis.h>
#include <math/optimize/efc.h>
#include <math/optimize/steep.h>
#include <math/optimize/update.h>
#include <math/optimize/mcsearch.h>

#include <math/scmat/linkage.h>

namespace sc {

static ForceLink<QNewtonOpt> optimize_force_link_a_;
static ForceLink<GDIISOpt> optimize_force_link_b_;
static ForceLink<EFCOpt> optimize_force_link_c_;
static ForceLink<BFGSUpdate> optimize_force_link_d_;
static ForceLink<PowellUpdate> optimize_force_link_e_;
static ForceLink<SteepestDescentOpt> optimize_force_link_f_;
static ForceLink<NewtonOpt> optimize_force_link_g_;
static ForceLink<MCSearch> optimize_force_link_h_;

}

#endif
