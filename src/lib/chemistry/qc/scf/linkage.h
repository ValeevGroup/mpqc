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

#ifndef _chemistry_qc_scf_linkage_h
#define _chemistry_qc_scf_linkage_h

#include <mpqc_config.h>

#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/scf/fbclhf.h>
#include <chemistry/qc/scf/hsoshf.h>
#include <chemistry/qc/scf/osshf.h>
#include <chemistry/qc/scf/tchf.h>
#include <chemistry/qc/scf/uhf.h>
#include <chemistry/qc/scf/scf.h>
#ifdef MPQC_NEW_FEATURES
#  include <chemistry/qc/scf/cadf/cadfclhf.h>
#  include <chemistry/qc/scf/cadf/approx_pairs.h>
#  include <chemistry/qc/scf/iter_logger.h>
#endif // MPQC_NEW_FEATURES

#include <math/scmat/linkage.h>
#include <chemistry/molecule/linkage.h>
#include <chemistry/qc/wfn/linkage.h>

#ifdef HAVE_LIBINT2
#  include <chemistry/qc/libint2/linkage.h>
#endif

#if MPQC_NEW_FEATURES && MPQC_HAS_ELEMENTAL
#  include<chemistry/qc/scf/taclhf.hpp>
#endif

namespace sc {

ForceLink<CLHF> scf_force_link_a_;
ForceLink<HSOSHF> scf_force_link_b_;
ForceLink<OSSHF> scf_force_link_c_;
ForceLink<TCHF> scf_force_link_d_;
ForceLink<UHF> scf_force_link_e_;
ForceLink<FockBuildCLHF> scf_force_link_f_;

#ifdef MPQC_NEW_FEATURES
ForceLink<ApproximatePairWriter> scf_force_link_g_;
ForceLink<SCFIterationLogger> scf_force_link_h_;
ForceLink<CADFCLHF> scf_force_link_i_;
#endif

#if MPQC_NEW_FEATURES && MPQC_HAS_ELEMENTAL
  ForceLink<mpqc::TA::CLHF> scf_force_link_g_;
#endif

} // namespace sc

#endif
