//
// localdf_runtime.h
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

#ifndef _mpqc_src_lib_chemistry_qc_lcao_localdf_runtime_h
#define _mpqc_src_lib_chemistry_qc_lcao_localdf_runtime_h

#include <chemistry/qc/lcao/df_runtime.h>
#include <chemistry/qc/lcao/tbint_runtime.h>

namespace sc {

  class LocalDensityFittingRuntime : public DensityFittingRuntimeBase {
    public:
      typedef LocalDensityFittingRuntime this_type;

      LocalDensityFittingRuntime(StateIn& si);

      void save_data_state(StateOut& so);

    protected:

      typedef Eigen::HouseholderQR<Eigen::MatrixXd> Decomposition;
      typedef std::map<std::pair<int, int>, std::shared_ptr<Decomposition> > DecompositionMap;
      DecompositionMap decomps_;


  };

} // end namespace sc

#endif
