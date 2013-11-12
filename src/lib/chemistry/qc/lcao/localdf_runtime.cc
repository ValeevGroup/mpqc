//
// localdf_runtime.cc
//
// Copyright (C) 2013 David Hollman
//
// Author: David Hollman <dhollman@vt.edu>
// Maintainer: DSH
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

#include <chemistry/qc/lcao/localdf_runtime.h>

using namespace sc;


static ClassDesc LocalDensityFittingRuntime_cd(
  typeid(LocalDensityFittingRuntime),
  "LocalDensityFittingRuntime",
  /* version = */ 1,
  "public DensityFittingRuntimeBase",
  /* pointer to default constructor = */ 0,
  /* pointer to keyval constructor = */ 0,
  /* pointer to StateIn constructor = */ create<LocalDensityFittingRuntime>
);


LocalDensityFittingRuntime::LocalDensityFittingRuntime(StateIn& si)
{
  throw FeatureNotImplemented("LocalDensityFittingRuntime construction from StateIn",
      __FILE__,
      __LINE__,
      class_desc()
  );
}

void
LocalDensityFittingRuntime::save_data_state(StateOut& so)
{
  throw FeatureNotImplemented("LocalDensityFittingRuntime writing to StateOut",
      __FILE__,
      __LINE__,
      class_desc()
  );

}
