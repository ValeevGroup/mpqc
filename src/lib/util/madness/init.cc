//
// init.cc
//
// Copyright (C) 2013 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <mpqc_config.h>

#ifdef HAVE_MADNESS
# ifndef WORLD_INSTANTIATE_STATIC_TEMPLATES
# define WORLD_INSTANTIATE_STATIC_TEMPLATES
# endif
# include <world/world.h>
#endif

#include <util/madness/init.h>
#include <util/misc/exenv.h>
#include <util/misc/scexception.h>

bool mpqc::MADNESSRuntime::mpqc_initialized_madness_ = false;

bool mpqc::MADNESSRuntime::initialized() {
  return madness::initialized();
}

void mpqc::MADNESSRuntime::initialize() {
  if (not madness::initialized()) {
    madness::initialize(sc::ExEnv::argc(), sc::ExEnv::argv());
    mpqc_initialized_madness_ = true;
  }
  else {
    mpqc_initialized_madness_ = false;
  }
  // now make sure that MADNESS has initialized MPI with full thread safety ...
  // if MPQC were using SafeMPI instead of MPI directly this would not be an issue
  if (SafeMPI::Query_thread() != MPI_THREAD_MULTIPLE && madness::World::get_default().rank() != 1) {
    throw sc::FeatureNotImplemented("nproc>1, and MPQC cannot get along with MADNESS because MADNESS was not configured with --with-mpi-thread=multiple; reconfigure MADNESS", __FILE__, __LINE__);
  }
}

void mpqc::MADNESSRuntime::finalize() {
  if (madness::initialized()) {
    if (mpqc_initialized_madness_)
      madness::finalize();
    mpqc_initialized_madness_ = false;
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
