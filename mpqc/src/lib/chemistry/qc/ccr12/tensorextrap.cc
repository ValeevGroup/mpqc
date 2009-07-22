#include <stdexcept>
#include <cassert>
#include <util/class/scexception.h>
#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/tensorextrap.h>


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


