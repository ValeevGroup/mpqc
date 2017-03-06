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

namespace detail {

template <typename T>
void print_cis_iteration(std::size_t iter, T norm, const EigenVector<T>& eig,
                         double time1, double time2) {
  ExEnv::out0() << indent << "iteration: " << iter << "\n";
  ExEnv::out0() << indent << "norm: " << norm << "\n";
  ExEnv::out0() << indent << "excitation energy: "
                << "\n";

  const auto size = eig.size();
  for (auto i = 0; i < size - 1; i++) {
    ExEnv::out0() << indent << indent << eig[i] << "\n";
  }
  ExEnv::out0() << indent << indent << eig[size - 1];

  ExEnv::out0() << "\n";
  ExEnv::out0() << indent << "total time: " << time1 + time2 << " S\n";
  ExEnv::out0() << indent << indent << "product time: " << time1 << " S\n";
  ExEnv::out0() << indent << indent << "davidson time: " << time2 << " S\n\n";
}
}

/**
 * CIS for closed shell system
 *
 */
template <typename Tile, typename Policy>
class CIS : public LCAOWavefunction<Tile, Policy>,
            public Provides<ExcitationEnergy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using numeric_type = typename Tile::numeric_type;

 private:
  /// precoditioner in DavidsonDiag, use F_aa - F_ii to
  /// approximate the diagonal of CIS H matrix
  struct Preconditioner {
    /// diagonal of F_ij matrix
    EigenVector<numeric_type> eps_o;
    /// diagonal of F_ab matrix
    EigenVector<numeric_type> eps_v;

    Preconditioner(const EigenVector<numeric_type>& eps_O,
                   const EigenVector<numeric_type>& eps_V)
        : eps_o(eps_O), eps_v(eps_V) {}

    void operator()(const numeric_type& e,
                    TA::DistArray<Tile, Policy>& guess) const {
      const auto& eps_v = this->eps_v;
      const auto& eps_o = this->eps_o;

      auto task = [&eps_v, &eps_o, e](TA::TensorD& result_tile) {
        const auto& range = result_tile.range();
        float norm = 0.0;
        for (const auto& i : range) {
          const auto result = result_tile[i] / (e + eps_o[i[0]] - eps_v[i[1]]);
          result_tile[i] = result;
          norm += result * result;
        }
        return std::sqrt(norm);
      };

      TA::foreach_inplace(guess, task);
      guess.world().gop.fence();
    }
  };

 public:
  /**
  * KeyVal constructor
  * @param kv
  * | Keyword | Type | Default| Description |
  * |---------|------|--------|-------------|
  * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
  * | max_iter| int | 30 | max number of iteration in davidson diagonalization|
  */
  explicit CIS(const KeyVal& kv) : LCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default Ref Wfn in CIS is not support! \n", __FILE__,
                       __LINE__, "ref");
    }
    max_iter_ = kv.value<int>("max_iter", 30);
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

  /// this approach computes and stores H matrix
  /// @return excitation energy
  std::vector<numeric_type> compute_cis(std::size_t n_roots,
                                        std::vector<TArray> guess_vector,
                                        double precision,
                                        bool triplets = false);

  /// @return guess vector of size n_roots as unit vector
  std::vector<TArray> init_guess_vector(std::size_t n_roots);

 private:
  std::size_t max_iter_;
  /// reference wavefunction
  std::shared_ptr<Wavefunction> ref_wfn_;
  /// eigen vector
  std::vector<TA::DistArray<Tile, Policy>> eigen_vector_;
  /// diagonal of f_ij
  EigenVector<numeric_type> eps_o_;
  /// diagonal of f_ab
  EigenVector<numeric_type> eps_v_;
};

template <typename Tile, typename Policy>
void CIS<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  if (!this->computed()) {
    auto& world = this->wfn_world()->world();
    auto target_precision = ex_energy->target_precision(0);
    auto target_ref_precision = target_precision / 100.0;

    this->init_sdref(ref_wfn_, target_ref_precision);

    auto time0 = mpqc::fenced_now(world);

    auto n_roots = ex_energy->n_roots();

    // get guess vector
    auto guess = init_guess_vector(n_roots);

    std::vector<numeric_type> result;

    if (ex_energy->singlets()) {
      result = compute_cis(n_roots, guess, target_precision);
    }

    // TODO separate singlets and triplets energy
    if (ex_energy->triplets()) {
      decltype(result) triplet_result =
          compute_cis(n_roots, guess, target_precision, true);
      result.insert(result.end(), triplet_result.begin(), triplet_result.end());
    }

    this->computed_ = true;
    this->set_value(ex_energy, result);

    auto time1 = mpqc::fenced_now(world);
    ExEnv::out0() << "CIS Total Time: " << mpqc::duration_in_s(time0, time1)
                  << " S\n";
  }
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::numeric_type>
CIS<Tile, Policy>::compute_cis(
    std::size_t n_roots, std::vector<typename CIS<Tile, Policy>::TArray> guess,
    double converge, bool triplets) {
  ExEnv::out0() << "\n";
  ExEnv::out0() << indent << "CIS: " << (triplets ? "Triplets" : "Singlets")
                << "\n";
  ExEnv::out0() << "\n";

  auto& world = this->wfn_world()->world();
  auto time0 = mpqc::fenced_now(world);

  auto& factory = this->lcao_factory();

  // compute required integrals
  auto F_ab = factory.compute(L"<a|F|b>");
  auto F_ij = factory.compute(L"<i|F|j>");
  auto I_ab = factory.compute(L"<a|I|b>");
  auto I_ij = factory.compute(L"<i|I|j>");
  auto G_iajb = factory.compute(L"<i a|G|j b>");

  // initialize diagonal
  if (eps_o_.size() == 0) {
    eps_o_ = array_ops::array_to_eigen(F_ij).diagonal();
  }
  if (eps_v_.size() == 0) {
    eps_v_ = array_ops::array_to_eigen(F_ab).diagonal();
  }

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

  auto time1 = mpqc::fenced_now(world);
  // used later
  auto time2 = mpqc::fenced_now(world);
  auto time = mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << indent << "Computed H matrix. Time: " << time << " S\n";

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, n_roots);

  auto pred = Preconditioner(eps_o_, eps_v_);

  // solve the lowest n_roots eigenvalues
  EigenVector<numeric_type> eig = EigenVector<numeric_type>::Zero(n_roots);
  auto i = 0;
  for (; i < max_iter_; i++) {
    time0 = mpqc::fenced_now(world);

    const auto n_v = guess.size();

    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    // product of H with guess vector
    for (auto i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      HB[i]("i,a") = H("i,j,a,b") * guess[i]("j,b");
    }

    time1 = mpqc::fenced_now(world);

    EigenVector<numeric_type> eig_new = dvd.extrapolate(HB, guess, pred);

    time2 = mpqc::fenced_now(world);

    auto norm = (eig - eig_new).norm();

    detail::print_cis_iteration(i, norm, eig_new,
                                mpqc::duration_in_s(time0, time1),
                                mpqc::duration_in_s(time1, time2));

    if (norm < converge) {
      break;
    }

    eig = eig_new;
  }

  ExEnv::out0() << "\n";

  if (i == max_iter_) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter_, "CIS");
  }

  eigen_vector_.insert(eigen_vector_.end(), guess.begin(), guess.end());

  return std::vector<numeric_type>(eig.data(), eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::TArray>
CIS<Tile, Policy>::init_guess_vector(std::size_t n_roots) {
  std::vector<TA::DistArray<Tile, Policy>> guess_vector(n_roots);

  auto& factory = this->lcao_factory();

  //  std::size_t n_i = factory.orbital_registry().retrieve("i").rank();
  std::size_t n_a = factory.orbital_registry().retrieve("a").rank();

  // use f_ia for shape
  auto f_ia = factory.compute(L"<i|F|a>");
  auto range = f_ia.trange();

  for (std::size_t i = 0; i < n_roots; i++) {
    TA::DistArray<Tile, Policy> guess(f_ia.world(), range, f_ia.shape());

    guess.fill(numeric_type(0.0));

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
      tile[element_idx] = numeric_type(1.0);
    }

    guess.truncate();

    guess_vector[i] = guess;
  }

  return guess_vector;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
