//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/integrals/lcao_factory.h>

namespace mpqc{
namespace qc{

class LCAOWfn : public Wfn {

  using ArrayType = typename Wfn::ArrayType;
  using LCAOFactoryType = integrals::LCAOFactory<TA::TensorD,TA::SparsePolicy>;

public:
  /*
   * KeyVal constructor
   * it includes all options from Wfn and LCAOFactory
   *
   */
  LCAOWfn(const KeyVal& kv);
  ~LCAOWfn() = default;

  LCAOFactoryType& lcao_factory();
  void compute(PropertyBase *pb) override;
  double value() override;

private:

  std::shared_ptr<LCAOFactoryType> lcao_factory_;

};

}
}


#endif //MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
