//
// Created by Chong Peng on 3/1/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_

#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/math/linalg/davidson_diag.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

/**
 * CIS for closed shell system
 *
 */
template <typename Tile, typename Policy>
class CIS : public LCAOWavefunction<Tile, Policy>,
            public Provides<ExcitationEnergy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  /**
  * KeyVal constructor
  * @param kv
  * | Keyword | Type | Default| Description |
  * |---------|------|--------|-------------|
  * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
  */
  explicit CIS(const KeyVal& kv) : LCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default Ref Wfn in CIS is not support! \n", __FILE__,
                       __LINE__, "ref");
    }
  }

  ~CIS() = default;

  std::vector<TArray> eigen_vector() const {
    if (eigen_vector_.empty()) {
      throw ProgrammingError("Eigenvector in CIS is not initialized!", __FILE__,
                             __LINE__);
    } else {
      return eigen_vector_;
    }
  }

 private:
  bool can_evaluate(ExcitationEnergy* ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy* ex_energy) override;

  /// @return excitation energy
  std::vector<typename Tile::numeric_type> compute_cis(std::size_t n_roots,
                                                       double converge,
                                                       bool triplets = false);

  /// @return guess vector of size n_roots as unit vector
  std::vector<TArray> init_guess_vector(std::size_t n_roots);

 private:
  /// reference wavefunction
  std::shared_ptr<Wavefunction> ref_wfn_;
  /// eigen vector
  std::vector<TA::DistArray<Tile, Policy>> eigen_vector_;
};

template <typename Tile, typename Policy>
void CIS<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  if (!this->computed()) {
    auto target_precision = ex_energy->target_precision(0);
    auto target_ref_precision = target_precision / 100.0;

    this->init_sdref(ref_wfn_, target_ref_precision);

    auto n_roots = ex_energy->n_roots();

    std::vector<typename Tile::numeric_type> result;

    if (ex_energy->singlets()) {
      result = compute_cis(n_roots, target_precision);
    }

    // TODO separate singlets and triplets energy
    if (ex_energy->triplets()) {
      decltype(result) triplet_result =
          compute_cis(n_roots, target_precision, true);
      result.insert(result.end(), triplet_result.begin(), triplet_result.end());
    }

    this->computed_ = true;
    this->set_value(ex_energy, result);
  }
}

template <typename Tile, typename Policy>
std::vector<typename Tile::numeric_type> CIS<Tile, Policy>::compute_cis(
    std::size_t n_roots, double converge, bool triplets) {
  auto& factory = this->lcao_factory();

  // compute required integrals
  auto F_ab = factory.compute(L"<a|F|b>");
  auto F_ij = factory.compute(L"<i|F|j>");
  auto I_ab = factory.compute(L"<a|I|b>");
  auto I_ij = factory.compute(L"<i|I|j>");
  auto G_iajb = factory.compute(L"<i a|G|j b>");

  // compute H
  TA::DistArray<Tile, Policy> H;

  // singlets
  if (!triplets) {
    auto G_ijab = factory.compute(L"<i j|G|a b>");
    H("i,j,a,b") = I_ij("i,j") * F_ab("a,b") - F_ij("i,j") * I_ab("a,b") +
                   2.0 * G_ijab("i,j,a,b") - G_iajb("i,a,j,b");
  }
  // triplets
  else {
    H("i,j,a,b") = I_ij("i,j") * F_ab("a,b") - F_ij("i,j") * I_ab("a,b") -
                   G_iajb("i,a,j,b");
  }

  factory.registry().purge_formula(L"<i j|G|a b>");
  factory.registry().purge_formula(L"<i a|G|j b>");

  // diagnoal elements of F_ab and F_ij
  EigenVector<typename Tile::numeric_type> eps_a, eps_i;
  {
    eps_a = array_ops::array_to_eigen(F_ab).diagonal();
    eps_i = array_ops::array_to_eigen(F_ij).diagonal();
  }

  // get guess vector
  auto guess = init_guess_vector(n_roots);

  // preconditioner
  // use f_aa - fii to estimate the diagonal of H
  auto pred = [&eps_a, &eps_i](const typename Tile::numeric_type& e,
                               TA::DistArray<Tile, Policy>& guess) {

    auto task = [&eps_a, &eps_i, &e](TA::TensorD& result_tile) {
      const auto& range = result_tile.range();
      float norm = 0.0;
      for (const auto& i : range) {
        const auto result = result_tile[i] / (e + eps_i[i[0]] - eps_a[i[1]]);
        result_tile[i] = result;
        norm += result * result;
      }
      return std::sqrt(norm);
    };

    TA::foreach_inplace(guess, task);
    guess.world().gop.fence();

  };

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, n_roots);

  const auto max_iter = 15;
  // solve the lowest n_roots eigenvalues
  EigenVector<typename Tile::numeric_type> eig =
      EigenVector<typename Tile::numeric_type>::Zero(n_roots);
  for (auto i = 0; i < max_iter; i++) {
    const auto n_v = guess.size();

    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    // product of H with guess vector
    for (auto i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      HB[i]("i,a") = H("i,j,a,b") * guess[i]("j,b");
    }

    EigenVector<typename Tile::numeric_type> eig_new =
        dvd.extrapolate(HB, guess, pred);

    auto norm = (eig - eig_new).norm();
    ExEnv::out0() << "Iter: " << i << "\n";
    ExEnv::out0() << "Norm: " << norm << "\n";
    ExEnv::out0() << "Energy: " << eig_new << "\n";

    if (norm < converge) {
      break;
    }

    eig = eig_new;
  }

  eigen_vector_.insert(eigen_vector_.end(), guess.begin(), guess.end());

  return std::vector<typename Tile::numeric_type>(eig.data(),
                                                  eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<TA::DistArray<Tile, Policy>> CIS<Tile, Policy>::init_guess_vector(
    std::size_t n_roots) {
  std::vector<TA::DistArray<Tile, Policy>> guess_vector(n_roots);

  auto& factory = this->lcao_factory();

  //  std::size_t n_i = factory.orbital_registry().retrieve("i").rank();
  std::size_t n_a = factory.orbital_registry().retrieve("a").rank();

  // use f_ia for shape
  auto f_ia = factory.compute(L"<i|F|a>");
  auto range = f_ia.trange();

  for (std::size_t i = 0; i < n_roots; i++) {
    TA::DistArray<Tile, Policy> guess(f_ia.world(), range, f_ia.shape());

    guess.fill(typename Tile::numeric_type(0.0));

    // fill in with 1.0
    std::size_t idx_i = i % n_a;
    std::size_t idx_a = i / n_a;

    //    std::cout << "element index" << std::endl;
    //    std::cout << idx_i << std::endl;
    //    std::cout << idx_a << std::endl;

    std::vector<std::size_t> element_idx{{idx_i, idx_a}};

    auto tile_idx = range.element_to_tile(element_idx);

    if (guess.is_local(tile_idx)) {
      // get tile
      auto tile = guess.find(tile_idx).get();
      // set value
      tile[element_idx] = typename Tile::numeric_type(1.0);
    }

    guess_vector[i] = guess;
  }

  return guess_vector;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
