#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

using namespace mpqc;

class MP2 : public lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>,
             public CanEvaluate<Energy> {
 public:
  MP2(const KeyVal& kv)
      : lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>(kv) {}
  virtual ~MP2() = default;

 protected:
  std::shared_ptr<lcao::Wavefunction> ref_wfn_;

  bool can_evaluate(Energy* energy) override { return energy->order() == 0; }

  void evaluate(Energy* energy) override { this->set_value(energy, 0); }
};

MPQC_CLASS_EXPORT2("MP2", MP2);
mpqc::detail::ForceLink<MP2> fl;
