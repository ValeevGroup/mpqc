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
#include <iostream>


//define this when GGspace = ip is to be used.
#define USE_GG_SPACE_EQ_IP

namespace sc {

class CCR12_Triples : virtual public RefCount {
  protected:
    CCR12_Info* z;

    Ref<Tensor> singles_intermediate_;
    Ref<Tensor> doubles_intermediate_;
    Ref<Tensor> rhs_intermediate_;
    Ref<Tensor> lhs_intermediate_;

    // Prediagonalization scheme
    // only for ii cases.
    void prediagon();
    void fill_in_ltensors();
    int pair_size_;

    RefDiagSCMatrix bdiag_;
    RefSCMatrix lmatrix_;
    Ref<Tensor> ltensor1_;
    Ref<Tensor> ltensor2_;

// Two cases. From here...
// for GGspace = ii
    void singles();
    void doubles();
    void denom_contraction();
    void denom_contraction_new();
    void offset_hhphhh(Ref<Tensor>&);
    double get_energy();
// for GGspace = ip
    void doubles_ig(Ref<Tensor>&, Ref<Tensor>&);
    void singles_ig(Ref<Tensor>&, Ref<Tensor>&) { std::cout << "== singles to be implemented ==" << std::endl;};
    void denom_contraction_ip();
    void offset_hgphhh(Ref<Tensor>&);
    double get_energy_ip();
// to here...

    void offset_bphhh(Ref<Tensor>&);

  public:
    CCR12_Triples(CCR12_Info* inz) : z(inz) {};
    ~CCR12_Triples() {};

    double compute();

};

}

#endif

