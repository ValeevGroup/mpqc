//
// ccr12_triples.h : computes unconventional triples correction to (T) model
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
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

#pragma once
#ifndef __chemistry_qc_ccr12_ccr12_triples_h
#define __chemistry_qc_ccr12_ccr12_triples_h

#include <chemistry/qc/ccr12/ccr12_info.h>

namespace sc {

class CCR12_Triples : virtual public RefCount {
  protected:
    CCR12_Info* z;

    Ref<Tensor> singles_intermediate_;
    Ref<Tensor> doubles_intermediate_;
    Ref<Tensor> rhs_intermediate_;
    Ref<Tensor> lhs_intermediate_;

    // for prediagonalization scheme
    void prediagon();
    void fill_in_ltensors();
    const int pair_size_;
    RefDiagSCMatrix bdiag_;
    RefSCMatrix lmatrix_;
    const RefSymmSCMatrix B_;
    const RefSymmSCMatrix X_;
    Ref<Tensor> ltensor1_;
    Ref<Tensor> ltensor2_;

    void singles();
    void doubles();
    void doubles_ig(Ref<Tensor>&, Ref<Tensor>&);
    void denom_contraction();
    void denom_contraction_new();
    void offset_hhphhh(Ref<Tensor>&);
    void offset_bphhh(Ref<Tensor>&);
    double get_energy();

  public:
    CCR12_Triples(CCR12_Info* inz, RefSymmSCMatrix b, RefSymmSCMatrix x)
     : z(inz), B_(b), X_(x), pair_size_(b.n()) {
      prediagon();
      singles_intermediate_ = new Tensor("singles_intermediate", z->mem());
      doubles_intermediate_ = new Tensor("doubles_intermediate", z->mem());
    };

    ~CCR12_Triples() {};

    double compute() {
      offset_hhphhh(singles_intermediate_);
      offset_hhphhh(doubles_intermediate_);
      singles(); // evaluating singles
      doubles(); // evaluating doubles
      singles_intermediate_->daxpy(doubles_intermediate_, 1.0); // adding doubles to singles to form lhs numerator
      rhs_intermediate_ = doubles_intermediate_->clone();
#if 0
      denom_contraction(); // contracting denominator to rhs numerator which is doubles
#else
      fill_in_ltensors();
      denom_contraction_new();
#endif

      return get_energy();
    };

};

}

#endif

