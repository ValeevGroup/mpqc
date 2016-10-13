//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

#include <string>

#include "../../../../../include/eigen.h"
#include "../../../../../utility/cc_utility.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/f12/cabs_singles.h>
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>
#include <mpqc/chemistry/qc/f12/f12_utility.h>
#include <mpqc/chemistry/qc/mbpt/mp2.h>

namespace mpqc {
namespace f12 {

class RMP2F12 : public qc::LCAOWavefunction{

public:
  using TArray = qc::LCAOWavefunction::ArrayType;
  using Matrix = RowMatrix<double>;

  RMP2F12(const KeyVal& kv);
  ~RMP2F12() = default;

  double value() override;
  std::tuple<Matrix,Matrix> compute();
  void compute(qc::PropertyBase* pb) override;
  void obsolete() override;

private:

  virtual TArray compute_B();
  virtual TArray compute_V();
  virtual TArray compute_X();
  virtual TArray compute_C();
  virtual std::tuple<TArray, TArray> compute_T();
  virtual double compute_cabs_singles();

protected:
  char approximation_;
  TA::SparseShape<float> ijij_ijji_shape_;

private:
  std::shared_ptr<qc::Wavefunction> ref_wfn_;
  double rmp2f12_energy_;
  bool cabs_singles_;
};

class RIRMP2F12 : public RMP2F12{

public:
  RIRMP2F12(const KeyVal& kv);
  ~RIRMP2F12() = default;

private:
  TArray compute_B() override ;
  TArray compute_V() override ;
  TArray compute_X() override ;
  TArray compute_C() override ;
  std::tuple<TArray, TArray> compute_T() override ;
  double compute_cabs_singles() override;

};

}  // end of namespace f12
}  // mpqc

#endif  // MPQC_MP2F12_H
