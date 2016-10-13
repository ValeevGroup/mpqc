//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DBMP2F12_H
#define MPQC_DBMP2F12_H

#include <mpqc/chemistry/qc/f12/db_f12_intermediates.h>
#include <mpqc/chemistry/qc/f12/mp2f12.h>
#include <mpqc/chemistry/qc/mbpt/dbmp2.h>

namespace mpqc {
namespace f12 {

class RIDBRMP2F12 : public qc::LCAOWavefunction{

public:

  using Matrix = RowMatrix<double>;

  RIDBRMP2F12(const KeyVal& kv);
  ~RIDBRMP2F12() = default;

  double value() override;
  double compute();
  void obsolete() override;
  void compute(qc::PropertyBase* pb) override;

private:
  std::tuple<Matrix,Matrix>  compute_db_mp2_f12_c();

private:
  std::shared_ptr<qc::Wavefunction> ref_wfn_; // scf
  double db_rmp2f12_energy_;

};
}  // namespace f12
}  // namespace mpqc

#endif  // MPQC_DBMP2F12_H
