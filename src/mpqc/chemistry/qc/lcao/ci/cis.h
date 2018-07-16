//
// Created by Chong Peng on 3/1/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_

#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/math/external/tiledarray/array_max_n.h"
#include "mpqc/math/linalg/davidson_diag.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/misc/print.h"

namespace mpqc {
namespace lcao {

namespace detail {

template <typename T>
struct IndexSort {
  std::pair<std::size_t, std::size_t> index;
  T value;

  IndexSort(std::size_t i, std::size_t j, T v)
      : index(std::make_pair(i, j)), value(v) {}

  IndexSort() = default;
  ~IndexSort() = default;

  bool operator<(const IndexSort &other) const { return value < other.value; }
};

}  // namespace detail

/**
 * CIS for closed shell system
 *
 * @warning This is not a efficient integral direct implementation of CIS, only
 * used to generate guess eigen vectors for EOM-CCSD
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
  struct Preconditioner : public DavidsonDiagPred<TArray> {
    /// diagonal of F_ij matrix
    EigenVector<numeric_type> eps_o;
    /// diagonal of F_ab matrix
    EigenVector<numeric_type> eps_v;

    Preconditioner(const EigenVector<numeric_type> &eps_O,
                   const EigenVector<numeric_type> &eps_V)
        : eps_o(eps_O), eps_v(eps_V) {}

    virtual void operator()(const EigenVector<numeric_type> &e,
                            std::vector<TArray> &guess) const {
      std::size_t n_roots = e.size();
      TA_ASSERT(n_roots == guess.size());
      for (std::size_t i = 0; i < n_roots; i++) {
        compute(e[i], guess[i]);
      }
    }

    void compute(const numeric_type &e,
                 TA::DistArray<Tile, Policy> &guess) const {
      const auto &eps_v = this->eps_v;
      const auto &eps_o = this->eps_o;

      auto task = [&eps_v, &eps_o, e](Tile &result_tile) {
        const auto &range = result_tile.range();
        float norm = 0.0;
        for (const auto &i : range) {
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
  // clang-format off
  /**
  * KeyVal constructor
  * @param kv
  * | Keyword | Type | Default| Description |
  * |---------|------|--------|-------------|
  * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
  * | method | df | standard or df | method to compute CIS, standard, df |
  * | max_iter| int | 30 | max number of iteration in davidson diagonalization|
  */
  // clang-format on
  explicit CIS(const KeyVal &kv) : LCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default Ref Wfn in CIS is not support! \n", __FILE__,
                       __LINE__, "ref");
    }
    max_iter_ = kv.value<int>("max_iter", 30);
    auto default_method =
        this->lcao_factory().basis_registry()->have(L"Κ") ? "df" : "standard";
    method_ = kv.value<std::string>("method", default_method);

    if (method_ != "df" && method_ != "standard" && method_ != "direct") {
      throw InputError("Invalid CIS method! \n", __FILE__, __LINE__, "method");
    }
    df_ = (method_ == "df" ? true : false);
  }

  ~CIS() = default;

  const std::vector<TArray> &eigen_vector() const {
    if (eigen_vector_.empty()) {
      throw ProgrammingError("Eigenvector in CIS is not initialized!", __FILE__,
                             __LINE__);
    } else {
      return eigen_vector_;
    }
  }

  const std::vector<numeric_type> eigen_value() const {
    if (eigen_value_.empty()) {
      throw ProgrammingError("Eigenvalue in CIS is not initialized!", __FILE__,
                             __LINE__);
    } else {
      return eigen_value_;
    }
  }

 private:
  bool can_evaluate(ExcitationEnergy *ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy *ex_energy) override;

  /// this approach stores two electron integral and computes and stores H
  /// matrix
  /// @return excitation energy
  std::vector<numeric_type> compute_cis(std::size_t n_roots,
                                        std::size_t n_guess, double precision,
                                        bool triplets = false);

  /// this approach uses density-fitting, it stores three center integral,
  /// but does not stores the H matrix, it compute the product of H with eigen
  /// vector
  /// @return excitation energy
  std::vector<numeric_type> compute_cis_df(std::size_t n_roots,
                                           std::size_t n_guess,
                                           double precision,
                                           bool triplets = false);

  /// @return guess vector of size n_roots as unit vector
  std::vector<TArray> init_guess_vector(std::size_t n_roots);

 private:
  /// max number of iteration in davidson
  std::size_t max_iter_;
  /// if has density fitting
  bool df_;
  /// CIS method string
  std::string method_;
  /// reference wavefunction
  std::shared_ptr<Wavefunction> ref_wfn_;
  /// eigen vector
  std::vector<TA::DistArray<Tile, Policy>> eigen_vector_;
  /// eigen value
  std::vector<numeric_type> eigen_value_;
  /// diagonal of f_ij
  EigenVector<numeric_type> eps_o_;
  /// diagonal of f_ab
  EigenVector<numeric_type> eps_v_;
};

template <typename Tile, typename Policy>
void CIS<Tile, Policy>::evaluate(ExcitationEnergy *ex_energy) {
  if (!this->computed()) {
    auto &world = this->wfn_world()->world();
    auto target_precision = ex_energy->target_precision(0);

    const auto &orbital_registry = this->ao_factory().orbital_registry();
    // if required orbitals not initialized, then initialize
    if (!orbital_registry.have(OrbitalIndex(L"i")) ||
        !orbital_registry.have(OrbitalIndex(L"a"))) {
      auto target_ref_precision = target_precision / 100.0;

      this->init_sdref(ref_wfn_, target_ref_precision);
    }

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nCIS Excitation Energy \n";
    auto n_roots = ex_energy->n_roots();
    auto n_guess = ex_energy->n_guess();

    std::vector<numeric_type> result;

    if (ex_energy->singlets()) {
      if (method_ == "standard") {
        result = compute_cis(n_roots, n_guess, target_precision);
      } else if (method_ == "df") {
        result = compute_cis_df(n_roots, n_guess, target_precision);
      }
    }

    // TODO separate singlets and triplets energy
    if (ex_energy->triplets()) {
      decltype(result) triplet_result;
      if (method_ == "standard") {
        triplet_result = compute_cis(n_roots, n_guess, target_precision, true);
      } else if (method_ == "df") {
        triplet_result =
            compute_cis_df(n_roots, n_guess, target_precision, true);
      }
      result.insert(result.end(), triplet_result.begin(), triplet_result.end());
    }

    this->computed_ = true;
    eigen_value_ = result;
    this->set_value(ex_energy, result);

    auto time1 = mpqc::fenced_now(world);
    ExEnv::out0() << "CIS Total Time: " << mpqc::duration_in_s(time0, time1)
                  << " S\n";
  }
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::numeric_type>
CIS<Tile, Policy>::compute_cis(std::size_t n_roots, std::size_t n_guess,
                               double converge, bool triplets) {
  ExEnv::out0() << "\n";
  ExEnv::out0() << indent
                << "CIS standard: " << (triplets ? "Triplets" : "Singlets")
                << "\n";
  ExEnv::out0() << "\n";

  auto &world = this->wfn_world()->world();
  auto time0 = mpqc::fenced_now(world);

  auto &factory = this->lcao_factory();

  // compute required integrals
  auto F_ab = factory.compute(L"<a|F|b>");
  auto F_ij = factory.compute(L"<i|F|j>");
  auto I_ab = factory.compute(L"<a|I|b>");
  auto I_ij = factory.compute(L"<i|I|j>");
  auto G_iajb = factory.compute(L"<i a|G|j b>");

  // initialize diagonal
  if (eps_o_.size() == 0) {
    eps_o_ = math::array_to_eigen(F_ij).diagonal();
  }
  if (eps_v_.size() == 0) {
    eps_v_ = math::array_to_eigen(F_ab).diagonal();
  }

  // get guess vector
  auto guess = init_guess_vector(n_guess);

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

  auto time1 = mpqc::fenced_now(world);
  auto time = mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << indent << "Computed H matrix. Time: " << time << " S\n";

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, true, 2, 10,
                                                10 * converge);

  auto pred = std::make_unique<Preconditioner>(eps_o_, eps_v_);

  auto oper = [&H](std::vector<TArray> &guess) {
    const auto n_v = guess.size();
    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    // product of H with guess vector
    for (std::size_t i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      HB[i]("i,a") = H("i,j,a,b") * guess[i]("j,b");
    }
    return HB;
  };
  // solve the lowest n_roots eigenvalues
  auto eig = dvd.solve(guess, oper, pred.get(), converge, max_iter_);

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, triplets);
  // get the latest eigen vector
  auto &eigen_vector = dvd.eigen_vector();

  for (std::size_t i = 0; i < n_roots; i++) {
    auto dominants = array_abs_max_n_index(eigen_vector[i], 5);

    ExEnv::out0() << "Dominant determinants of excited wave function " << i + 1
                  << "\n";
    util::print_cis_dominant_elements(dominants);
    ExEnv::out0() << "\n";
  }

  eigen_vector_.insert(eigen_vector_.end(), eigen_vector.begin(),
                       eigen_vector.end());

  return std::vector<numeric_type>(eig.data(), eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::numeric_type>
CIS<Tile, Policy>::compute_cis_df(std::size_t n_roots, std::size_t n_guess,
                                  double converge, bool triplets) {
  ExEnv::out0() << "\n";
  ExEnv::out0() << indent << "CIS Density-fitting: "
                << (triplets ? "Triplets" : "Singlets") << "\n";
  ExEnv::out0() << "\n";

  auto &factory = this->lcao_factory();

  // compute required integrals
  auto F_ab = factory.compute(L"<a|F|b>[df]");
  auto F_ij = factory.compute(L"<i|F|j>[df]");
  auto I_ab = factory.compute(L"<a|I|b>");
  auto I_ij = factory.compute(L"<i|I|j>");
  auto X_ab = factory.compute(L"(Κ|G|a b)[inv_sqr]");
  auto X_ij = factory.compute(L"(Κ|G|i j)[inv_sqr]");
  auto X_ia = factory.compute(L"(Κ|G|i a)[inv_sqr]");

  // initialize diagonal
  if (eps_o_.size() == 0) {
    eps_o_ = math::array_to_eigen(F_ij).diagonal();
  }
  if (eps_v_.size() == 0) {
    eps_v_ = math::array_to_eigen(F_ab).diagonal();
  }

  // get guess vector
  auto guess = init_guess_vector(n_guess);

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, true, 2, 10,
                                                10 * converge);

  auto pred = std::make_unique<Preconditioner>(eps_o_, eps_v_);

  auto oper = [&F_ab, &F_ij, &I_ab, &I_ij, &X_ab, &X_ij, &X_ia,
               triplets](std::vector<TArray> &guess) {
    const auto n_v = guess.size();
    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    for (std::size_t i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      const auto &vec = guess[i];
      // singlets
      if (!triplets) {
        HB[i]("j,b") = vec("i,a") * I_ij("i,j") * F_ab("a,b") -
                       F_ij("i,j") * vec("i,a") * I_ab("a,b") +
                       2.0 * vec("i,a") * X_ia("x,i,a") * X_ia("x,j,b") -
                       X_ab("x,a,b") * vec("i,a") * X_ij("x,i,j");
      }
      // triplets
      else {
        HB[i]("j,b") = vec("i,a") * I_ij("i,j") * F_ab("a,b") -
                       F_ij("i,j") * vec("i,a") * I_ab("a,b") +
                       -X_ab("x,a,b") * vec("i,a") * X_ij("x,i,j");
      }
    }
    return HB;
  };
  // solve the lowest n_roots eigenvalues
  auto eig = dvd.solve(guess, oper, pred.get(), converge, max_iter_);

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, triplets);
  // get the latest eigen vector
  auto &eigen_vector = dvd.eigen_vector();

  for (std::size_t i = 0; i < n_roots; i++) {
    auto dominants = array_abs_max_n_index(eigen_vector[i], 5);

    ExEnv::out0() << "Dominant determinants of excited wave function " << i + 1
                  << "\n";
    util::print_cis_dominant_elements(dominants);
    ExEnv::out0() << "\n";
  }

  eigen_vector_.insert(eigen_vector_.end(), eigen_vector.begin(),
                       eigen_vector.end());

  return std::vector<numeric_type>(eig.data(), eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::TArray>
CIS<Tile, Policy>::init_guess_vector(std::size_t n_roots) {
  std::vector<TA::DistArray<Tile, Policy>> guess_vector(n_roots);

  auto &factory = this->lcao_factory();
  auto& world = factory.world();

  std::size_t n_m = factory.orbital_registry().retrieve("m").rank();
  const auto &orbital_i = factory.orbital_registry().retrieve("i");
  const auto &orbital_a = factory.orbital_registry().retrieve("a");
  std::size_t n_i = orbital_i.rank();
  std::size_t n_a = orbital_a.rank();
  std::size_t n_frozen = n_m - n_i;

  auto canonical_orbs_provider = std::dynamic_pointer_cast<
      typename CanonicalOrbitalSpace<TArray>::Provider>(ref_wfn_);
  bool has_canonical_orbs =
      canonical_orbs_provider && canonical_orbs_provider->can_evaluate();

  // these two will be computed if is not canonical
  TArray local_to_canonical_occ;
  TArray local_to_canonical_unocc;

  EigenVector<numeric_type> energy_occ;
  EigenVector<numeric_type> energy_unocc;
  if (has_canonical_orbs) {
    energy_occ = eps_o_;
    energy_unocc = eps_v_;
  } else {
    auto f_pq =
        df_ ? factory.compute(L"<p|F|q>[df]") : factory.compute(L"<p|F|q>");

    auto f_eig = math::array_to_eigen(f_pq);
    Eigen::VectorXd evals;
    RowMatrixXd C;
    if (world.rank() == 0) {
      Eigen::SelfAdjointEigenSolver<decltype(f_eig)> es(f_eig);
      evals = es.eigenvalues();
      C = es.eigenvectors();
    }
    world.gop.broadcast_serializable(evals, 0);
    world.gop.broadcast_serializable(C, 0);

    energy_occ = evals.segment(n_frozen, n_i);
    energy_unocc = evals.segment(n_m, n_a);

    decltype(C) C_i = C.block(n_frozen, n_frozen, n_i, n_i);
    decltype(C) C_a = C.block(n_m, n_m, n_a, n_a);

    auto &tr_i = orbital_i.trange();
    auto &tr_a = orbital_a.trange();
    auto &world = f_pq.world();

    local_to_canonical_occ =
        math::eigen_to_array<Tile, Policy>(world, C_i, tr_i, tr_i);
    local_to_canonical_unocc =
        math::eigen_to_array<Tile, Policy>(world, C_a, tr_a, tr_a);
  }
  // use f_ia for shape
  auto f_ia =
      df_ ? factory.compute(L"<i|F|a>[df]") : factory.compute(L"<i|F|a>");
  auto range = f_ia.trange();

  // find the minimum n_roots of F_ii - F_aa
  std::vector<detail::IndexSort<numeric_type>> index(n_i * n_a);
  {
    // TODO might want to parallel this section
    for (std::size_t i = 0; i < n_i; ++i) {
      for (std::size_t a = 0; a < n_a; ++a) {
        index[i * n_a + a] = std::move(detail::IndexSort<numeric_type>(
            i, a, energy_unocc[a] - energy_occ[i]));
      }
    }

    std::sort(index.begin(), index.end());

    //    for (std::size_t i = 0; i < n_roots; i++) {
    //      std::cout << "index: " << index[i].index.first << " "
    //                << index[i].index.second << " " << index[i].value << "\n";
    //    }
  }

  for (std::size_t i = 0; i < n_roots; i++) {
    TA::DistArray<Tile, Policy> guess(f_ia.world(), range, f_ia.shape());

    guess.fill(numeric_type(0.0));

    std::size_t idx_i = index[i].index.first;
    std::size_t idx_a = index[i].index.second;
    //    std::size_t idx_i = n_i - 1;
    //    std::size_t idx_a = i;
    //        std::cout << "element index" << std::endl;
    //        std::cout << idx_i << std::endl;
    //        std::cout << idx_a << std::endl;

    // fill in with 1.0
    std::vector<std::size_t> element_idx{{idx_i, idx_a}};

    auto tile_idx = range.element_to_tile(element_idx);

    if (guess.is_local(tile_idx)) {
      // get tile
      auto tile = guess.find(tile_idx).get();
      // set value
      tile[element_idx] = numeric_type(1.0);
    }

    guess.truncate();
    // convert to local basis
    if (!has_canonical_orbs) {
      guess("i,a") = guess("i1, a1") * local_to_canonical_occ("i,i1") *
                     local_to_canonical_unocc("a,a1");
    }

    guess_vector[i] = guess;
  }

  return guess_vector;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
