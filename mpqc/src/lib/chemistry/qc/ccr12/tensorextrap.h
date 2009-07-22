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
