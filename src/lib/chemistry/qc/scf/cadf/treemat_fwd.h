//
// treemat_fwd.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: May 2, 2014
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

#ifndef _chemistry_qc_scf_cadf_treemat_fwd_h
#define _chemistry_qc_scf_cadf_treemat_fwd_h

#include <Eigen/Dense>

namespace sc { namespace cadf {

typedef uint64_t uli;
typedef unsigned int uint;

template<
  typename NormContainer=Eigen::VectorXd,
  typename Index=uli,
  typename NormValue=typename Eigen::internal::traits<NormContainer>::Scalar
>
class TreeBlock;

template<typename BlockType=TreeBlock<>>
class TreeMatrix;


}} // end namespaces

#endif /* _chemistry_qc_scf_cadf_treemat_fwd_h */
