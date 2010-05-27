//
// ccr12_triples.cc
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@theochem.uni-stuttgart.de>
// Maintainer: TS
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


#include <chemistry/qc/ccr12/ccr12_triples.h>
#include <iostream>

using namespace sc;
using namespace std;

double CCR12_Triples::compute() {


#ifdef USE_GG_SPACE_EQ_IP

  singles_intermediate_ = new Tensor("singles_intermediate", z->mem());
  offset_hgphhh(singles_intermediate_);
  doubles_intermediate_ = singles_intermediate_->clone();
  rhs_intermediate_ = doubles_intermediate_->clone();

  // evaluating V * t2 and V * t1
  doubles_ig(doubles_intermediate_);
  singles_ig(singles_intermediate_);
  singles_intermediate_->daxpy(doubles_intermediate_, 1.0); // adding doubles to singles to form lhs numerator

  rhs_intermediate_ = doubles_intermediate_->clone();
  denom_contraction_ig(); // contracting denominator to rhs numerator which is doubles

  return get_energy_ig();

#else

  std::cout << "== Using obsolete iip/iii ansatz ==" << std::endl;
  prediagon();
  singles_intermediate_ = new Tensor("singles_intermediate", z->mem());
  doubles_intermediate_ = new Tensor("doubles_intermediate", z->mem());
  offset_hhphhh(singles_intermediate_);
  offset_hhphhh(doubles_intermediate_);
  singles(); // evaluating singles
  doubles(); // evaluating doubles
  singles_intermediate_->daxpy(doubles_intermediate_, 1.0); // adding doubles to singles to form lhs numerator
  rhs_intermediate_ = doubles_intermediate_->clone();

// this toggles whether we use prediagonalization scheme or not.
#if 0
  denom_contraction(); // contracting denominator to rhs numerator which is doubles
#else
  fill_in_ltensors();
  denom_contraction_new();
#endif

  return get_energy();

#endif
};
