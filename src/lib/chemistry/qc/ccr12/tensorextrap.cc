//
// tensorextrap.cc
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

#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/tensorextrap.h>
#include <util/class/scexception.h>


using namespace sc;
using namespace std;

static ClassDesc TensorExtrapData_cd(
  typeid(TensorExtrapData),"TensorExtrapData",1,"public SCExtrapData",
  0, 0, create<TensorExtrapData>);

TensorExtrapData::TensorExtrapData(StateIn& s): SCExtrapData(s){
  throw ProgrammingError("TensorExtrapData::TensorExtrapData(StateIn& s) not implemented",__FILE__,__LINE__);
}


TensorExtrapData::TensorExtrapData(const Ref<Tensor>& mat){
  m=mat;
}


void TensorExtrapData::save_data_state(StateOut& s){
  throw ProgrammingError("TensorExtrapData::save_data_state(StateOut& s) not implemented",__FILE__,__LINE__);
}


SCExtrapData* TensorExtrapData::copy(){
  return new TensorExtrapData(m->copy());
}


void TensorExtrapData::zero(){
  m->zero();
}


void TensorExtrapData::accumulate_scaled(double scale,const Ref<SCExtrapData>& data){

  TensorExtrapData* a=require_dynamic_cast<TensorExtrapData*>(data.pointer(),
                                      "TensorExtrapData::accumulate_scaled");
  m->daxpy(a->m,scale); 
}


/////////////////////////////////////////////////////////////////////////////////////

static ClassDesc TensorExtrapError_cd(
  typeid(TensorExtrapError),"TensorExtrapError",1,"public SCExtrapError",
  0, 0, create<TensorExtrapError>);


TensorExtrapError::TensorExtrapError(StateIn& s): SCExtrapError(s){
  throw ProgrammingError("TensorExtrapError::TensorExtrapError(StateIn& s) not implemented",__FILE__,__LINE__);
}


TensorExtrapError::TensorExtrapError(const Ref<Tensor>& mat){
  m=mat;
}


void TensorExtrapError::save_data_state(StateOut& s){
  throw ProgrammingError("TensorExtrapError::save_data_state(StateOut& s) not implemented",__FILE__,__LINE__);
}


double TensorExtrapError::error(){
  return 0.0; // perhaps it is not needed  
}


double TensorExtrapError::scalar_product(const Ref<SCExtrapError>& arg){
  TensorExtrapError* a=require_dynamic_cast<TensorExtrapError*>(arg.pointer(),"TensorExtrapError::scalar_product");
  return m->ddot(a->m);
} 


