//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/integrals/lcao_factory.h>
#include "../../../../../utility/trange1_engine.h"

namespace mpqc{
namespace qc{

class LCAOWavefunction : public Wavefunction {

  using ArrayType = typename Wavefunction::ArrayType;
  using LCAOFactoryType = integrals::LCAOFactory<TA::TensorD,TA::SparsePolicy>;

public:
  /*
   * KeyVal constructor
   * it includes all options from Wfn and LCAOFactory
   *
   */
  LCAOWavefunction(const KeyVal& kv);
  ~LCAOWavefunction() = default;

  LCAOFactoryType& lcao_factory();
  void compute(PropertyBase *pb) override;
  double value() override;
  void obsolete() override;

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return trange1_engine_;
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return orbital_energy_;
  }

  bool is_frozen_core() const;
  size_t occ_block() const;
  size_t unocc_block() const;

protected:
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;
  std::shared_ptr<mpqc::TRange1Engine> trange1_engine_;

private:

  std::shared_ptr<LCAOFactoryType> lcao_factory_;
  bool frozen_core_;
  std::size_t occ_block_;
  std::size_t unocc_block_;

};

}
}


#endif //MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
