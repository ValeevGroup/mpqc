//
// tensorextrap.h --- interface to SCExtrapData and SCExtrapError 
//
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

#ifndef _chemistry_qc_ccr12_tensorextrap_h
#define _chemistry_qc_ccr12_tensorextrap_h

#include <math/optimize/scextrap.h>
#include <chemistry/qc/ccr12/tensor.h>

namespace sc {

class TensorExtrapData: public SCExtrapData {
  private:
    Ref<Tensor> m;
  public:
    TensorExtrapData(StateIn&);
    TensorExtrapData(const Ref<Tensor>&);

    void save_data_state(StateOut&);

    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const Ref<SCExtrapData>&);
};


class TensorExtrapError: public SCExtrapError {
  private:
    Ref<Tensor> m;
  public:
    TensorExtrapError(StateIn&);
    TensorExtrapError(const Ref<Tensor>&);

    void save_data_state(StateOut&);

    double error();    
    double scalar_product(const Ref<SCExtrapError>&); 

};

}
#endif

