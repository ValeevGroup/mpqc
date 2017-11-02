//
// Created by Chong Peng on 7/21/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_IP_EOM_CCSD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_IP_EOM_CCSD_H_

#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/math/linalg/davidson_diag.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class IP_EOM_CCSD : public CCSD<Tile, Policy>,
                    public Provides<ExcitationEnergy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using GuessVector = ::mpqc::cc::T1T2<TArray, TArray>;
  using numeric_type = typename Tile::numeric_type;

  // clang-format off
  /**
  * KeyVal constructor
  * @param kv
  *
  * | Keyword | Type | Default| Description |
  * |---------|------|--------|-------------|
  * | max_vector | int | 8 | max number of guess vector per root |
  * | vector_threshold | double | 1.0e-5 | threshold for the norm of new guess vector |
  */

  // clang-format on

  IP_EOM_CCSD(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    max_vector_ = kv.value<int>("max_vector", 8);
    vector_threshold_ = kv.value<double>("vector_threshold", 1.0e-5);
  }

  void obsolete() override { CCSD<Tile, Policy>::obsolete(); }

 protected:
  using CCSD<Tile, Policy>::can_evaluate;
  using CCSD<Tile, Policy>::evaluate;

  bool can_evaluate(ExcitationEnergy *ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy *ex_energy) override;

 private:
  /// @return guess vector of size n_roots as unit vector
  std::vector<GuessVector> init_guess_vector(std::size_t n_roots);

  /// @return ionization potentials
  EigenVector<numeric_type> ip_eom_ccsd_davidson_solver(
      std::vector<GuessVector> &C, const cc::Intermediates<TArray> &imds,
      std::size_t max_iter, double convergence);

  /// compute the product of H with vectors
  TArray compute_HS1(const TArray &Ci, const TArray &Caij,
                     const cc::Intermediates<TArray> &imds);
  TArray compute_HS2(const TArray &Ci, const TArray &Caij,
                     const cc::Intermediates<TArray> &imds);

 private:
  std::size_t max_vector_;   // max number of guess vector
  double vector_threshold_;  // threshold for norm of new guess vector
};

#if TA_DEFAULT_POLICY == 0
extern template class IP_EOM_CCSD<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class IP_EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
#endif

}  // end of namespace lcao
}  // end of namespace mpqc

#include "ip_eom_ccsd_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_IP_EOM_CCSD_H_
