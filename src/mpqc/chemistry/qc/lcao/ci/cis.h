//
// Created by Chong Peng on 3/1/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CI_CIS_H_

#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/math/linalg/davidson_diag.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/misc/print.h"

namespace mpqc {
namespace lcao {

namespace detail {}  // namespace detail

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
  struct Preconditioner : public DavidsonDiagPreconditioner<TArray> {
    /// diagonal of F_ij matrix
    EigenVector<numeric_type> eps_o;
    /// diagonal of F_ab matrix
    EigenVector<numeric_type> eps_v;

    Preconditioner(const EigenVector<numeric_type> &eps_O,
                   const EigenVector<numeric_type> &eps_V)
        : eps_o(eps_O), eps_v(eps_V) {}

    using DavidsonDiagPreconditioner<TArray>::operator();

    virtual void compute(const numeric_type &e,
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
  * | method | df | standard or df | method to compute CIS, standard, df or direct |
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

  std::vector<TArray> eigen_vector() const {
    if (eigen_vector_.empty()) {
      throw ProgrammingError("Eigenvector in CIS is not initialized!", __FILE__,
                             __LINE__);
    } else {
      return eigen_vector_;
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
                                        std::vector<TArray> guess_vector,
                                        double precision,
                                        bool triplets = false);

  /// this approach uses density-fitting, it stores three center integral,
  /// but does not stores the H matrix, it compute the product of H with eigen
  /// vector
  /// @return excitation energy
  std::vector<numeric_type> compute_cis_df(std::size_t n_roots,
                                           std::vector<TArray> guess_vector,
                                           double precision,
                                           bool triplets = false);

  /// this approach is integral direct, it does not stores the H matrix and two
  /// electron integral, it compute the product of H with eigen vector using
  /// direct AO integral
  /// @return excitation energy
  std::vector<numeric_type> compute_cis_direct(std::size_t n_roots,
                                               std::vector<TArray> guess_vector,
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
    auto target_ref_precision = target_precision / 100.0;

    this->init_sdref(ref_wfn_, target_ref_precision);

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nCIS Excitation Energy \n";
    auto n_roots = ex_energy->n_roots();

    // get guess vector
    auto guess = init_guess_vector(n_roots);

    std::vector<numeric_type> result;

    if (ex_energy->singlets()) {
      if (method_ == "standard") {
        result = compute_cis(n_roots, guess, target_precision);
      } else if (method_ == "df") {
        result = compute_cis_df(n_roots, guess, target_precision);
      } else if (method_ == "direct") {
        result = compute_cis_direct(n_roots, guess, target_precision);
      }
    }

    // TODO separate singlets and triplets energy
    if (ex_energy->triplets()) {
      decltype(result) triplet_result;
      if (method_ == "standard") {
        triplet_result = compute_cis(n_roots, guess, target_precision, true);
      } else if (method_ == "df") {
        triplet_result = compute_cis_df(n_roots, guess, target_precision, true);
      } else if (method_ == "direct") {
        triplet_result =
            compute_cis_direct(n_roots, guess, target_precision, true);
      }
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

  auto time1 = mpqc::fenced_now(world);
  // used later
  auto time2 = mpqc::fenced_now(world);
  auto time = mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << indent << "Computed H matrix. Time: " << time << " S\n";

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, true, 2, 10);

  auto pred = Preconditioner(eps_o_, eps_v_);

  // solve the lowest n_roots eigenvalues
  EigenVector<numeric_type> eig = EigenVector<numeric_type>::Zero(n_roots);
  std::size_t i = 0;
  for (; i < max_iter_; i++) {
    time0 = mpqc::fenced_now(world);

    const auto n_v = guess.size();

    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    // product of H with guess vector
    for (std::size_t i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      HB[i]("i,a") = H("i,j,a,b") * guess[i]("j,b");
    }

    time1 = mpqc::fenced_now(world);

    EigenVector<numeric_type> eig_new = dvd.extrapolate(HB, guess, pred);

    time2 = mpqc::fenced_now(world);

    EigenVector<numeric_type> delta_e = eig - eig_new;
    auto norm = delta_e.norm();

    util::print_excitation_energy_iteration(i, delta_e, eig_new,
                                            mpqc::duration_in_s(time0, time1),
                                            mpqc::duration_in_s(time1, time2));

    if (norm < converge) {
      break;
    }

    eig = eig_new;
  }

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, triplets);

  if (i == max_iter_) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter_, "CIS");
  }

  // get the latest eigen vector
  auto &eigen_vector = dvd.eigen_vector().back();
  eigen_vector_.insert(eigen_vector_.end(), eigen_vector.begin(),
                       eigen_vector.end());

  return std::vector<numeric_type>(eig.data(), eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::numeric_type>
CIS<Tile, Policy>::compute_cis_df(
    std::size_t n_roots, std::vector<typename CIS<Tile, Policy>::TArray> guess,
    double converge, bool triplets) {
  ExEnv::out0() << "\n";
  ExEnv::out0() << indent << "CIS Density-fitting: "
                << (triplets ? "Triplets" : "Singlets") << "\n";
  ExEnv::out0() << "\n";

  auto &world = this->wfn_world()->world();
  auto &factory = this->lcao_factory();
  auto &ao_factory = this->ao_factory();

  // compute required integrals
  auto F_ab = factory.compute(L"<a|F|b>[df]");
  auto F_ij = factory.compute(L"<i|F|j>[df]");
  auto I_ab = factory.compute(L"<a|I|b>");
  auto I_ij = factory.compute(L"<i|I|j>");
  auto X_ab = factory.compute(L"(Κ|G|a b)");
  auto X_ij = factory.compute(L"(Κ|G|i j)");
  auto X_ia = factory.compute(L"(Κ|G|i a)");
  auto X = ao_factory.compute(L"(Κ|G|Λ)[inv]");

  // initialize diagonal
  if (eps_o_.size() == 0) {
    eps_o_ = array_ops::array_to_eigen(F_ij).diagonal();
  }
  if (eps_v_.size() == 0) {
    eps_v_ = array_ops::array_to_eigen(F_ab).diagonal();
  }

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, true, 2, 10);

  auto pred = Preconditioner(eps_o_, eps_v_);

  // solve the lowest n_roots eigenvalues
  EigenVector<numeric_type> eig = EigenVector<numeric_type>::Zero(n_roots);
  std::size_t i = 0;
  for (; i < max_iter_; i++) {
    auto time0 = mpqc::fenced_now(world);

    const auto n_v = guess.size();

    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    // product of H with guess vector
    for (std::size_t i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      const auto &vec = guess[i];
      // singlets
      if (!triplets) {
        HB[i]("j,b") =
            vec("i,a") * I_ij("i,j") * F_ab("a,b") -
            F_ij("i,j") * vec("i,a") * I_ab("a,b") +
            2.0 * vec("i,a") * X_ia("x,i,a") * X("x,y") * X_ia("y,j,b") -
            X_ab("x,a,b") * vec("i,a") * X("x,y") * X_ij("y,i,j");
      }
      // triplets
      else {
        HB[i]("j,b") = vec("i,a") * I_ij("i,j") * F_ab("a,b") -
                       F_ij("i,j") * vec("i,a") * I_ab("a,b") +
                       -X_ab("x,a,b") * vec("i,a") * X("x,y") * X_ij("y,i,j");
      }
    }

    auto time1 = mpqc::fenced_now(world);

    EigenVector<numeric_type> eig_new = dvd.extrapolate(HB, guess, pred);

    auto time2 = mpqc::fenced_now(world);

    EigenVector<numeric_type> delta_e = eig - eig_new;
    auto norm = delta_e.norm();

    util::print_excitation_energy_iteration(i, delta_e, eig_new,
                                            mpqc::duration_in_s(time0, time1),
                                            mpqc::duration_in_s(time1, time2));

    if (norm < converge) {
      break;
    }

    eig = eig_new;
  }

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, triplets);

  if (i == max_iter_) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter_, "CIS");
  }

  // get the latest eigen vector
  auto &eigen_vector = dvd.eigen_vector().back();

  eigen_vector_.insert(eigen_vector_.end(), eigen_vector.begin(),
                       eigen_vector.end());

  return std::vector<numeric_type>(eig.data(), eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::numeric_type>
CIS<Tile, Policy>::compute_cis_direct(
    std::size_t n_roots, std::vector<typename CIS<Tile, Policy>::TArray> guess,
    double converge, bool triplets) {
  ExEnv::out0() << "\n";
  ExEnv::out0() << indent
                << "CIS Direct: " << (triplets ? "Triplets" : "Singlets")
                << "\n";
  ExEnv::out0() << "\n";

  auto &world = this->wfn_world()->world();
  auto &factory = this->lcao_factory();
  auto &ao_factory = this->ao_factory();

  // compute required integrals
  auto F_ab = factory.compute(L"<a|F|b>");
  auto F_ij = factory.compute(L"<i|F|j>");
  auto I_ab = factory.compute(L"<a|I|b>");
  auto I_ij = factory.compute(L"<i|I|j>");
  auto G = ao_factory.compute_direct(L"(μ ν| G|κ λ)");
  auto C_a = factory.orbital_registry().retrieve(L"a");
  auto C_i = factory.orbital_registry().retrieve(L"i");

  // initialize diagonal
  if (eps_o_.size() == 0) {
    eps_o_ = array_ops::array_to_eigen(F_ij).diagonal();
  }
  if (eps_v_.size() == 0) {
    eps_v_ = array_ops::array_to_eigen(F_ab).diagonal();
  }

  // davidson object
  DavidsonDiag<TA::DistArray<Tile, Policy>> dvd(n_roots, true, 2, 10);

  auto pred = Preconditioner(eps_o_, eps_v_);

  // solve the lowest n_roots eigenvalues
  EigenVector<numeric_type> eig = EigenVector<numeric_type>::Zero(n_roots);
  std::size_t i = 0;
  for (; i < max_iter_; i++) {
    auto time0 = mpqc::fenced_now(world);

    const auto n_v = guess.size();

    std::vector<TA::DistArray<Tile, Policy>> HB(n_v);
    // TODO need to avoid compute all four center at each iteration
    // product of H with guess vector
    for (std::size_t i = 0; i < n_v; i++) {
      //    std::cout << guess[i] << std::endl;
      const auto &vec = guess[i];
      // singlets
      if (!triplets) {
        HB[i]("j,b") =
            vec("i,a") * I_ij("i,j") * F_ab("a,b") -
            F_ij("i,j") * vec("i,a") * I_ab("a,b") +
            2.0 *
                ((vec("i,a") * C_a("nu,a") * C_i("mu,i")) *
                 G("mu,nu,rho,sigma")) *
                C_i("rho,j") * C_a("sigma,b") -
            ((vec("i,a") * C_a("rho,a") * C_i("mu,i")) * G("mu,nu,rho,sigma")) *
                C_i("nu,j") * C_a("sigma,b");
      }
      // triplets
      else {
        HB[i]("j,b") = vec("i,a") * I_ij("i,j") * F_ab("a,b") -
                       F_ij("i,j") * vec("i,a") * I_ab("a,b") +
                       -((vec("i,a") * C_a("rho,a") * C_i("mu,i")) *
                         G("mu,nu,rho,sigma")) *
                           C_i("nu,j") * C_a("sigma,b");
      }
    }

    auto time1 = mpqc::fenced_now(world);

    EigenVector<numeric_type> eig_new = dvd.extrapolate(HB, guess, pred);

    auto time2 = mpqc::fenced_now(world);

    EigenVector<numeric_type> delta_e = eig - eig_new;
    auto norm = delta_e.norm();

    util::print_excitation_energy_iteration(i, delta_e, eig_new,
                                            mpqc::duration_in_s(time0, time1),
                                            mpqc::duration_in_s(time1, time2));

    if (norm < converge) {
      break;
    }

    eig = eig_new;
  }

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, triplets);

  if (i == max_iter_) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter_, "CIS");
  }

  // get the latest eigen vector
  auto &eigen_vector = dvd.eigen_vector().back();

  eigen_vector_.insert(eigen_vector_.end(), eigen_vector.begin(),
                       eigen_vector.end());

  return std::vector<numeric_type>(eig.data(), eig.data() + eig.size());
}

template <typename Tile, typename Policy>
std::vector<typename CIS<Tile, Policy>::TArray>
CIS<Tile, Policy>::init_guess_vector(std::size_t n_roots) {
  std::vector<TA::DistArray<Tile, Policy>> guess_vector(n_roots);

  auto &factory = this->lcao_factory();

  std::size_t n_i = factory.orbital_registry().retrieve("i").rank();
  //  std::size_t n_a = factory.orbital_registry().retrieve("a").rank();

  // use f_ia for shape
  auto f_ia =
      df_ ? factory.compute(L"<i|F|a>[df]") : factory.compute(L"<i|F|a>");
  auto range = f_ia.trange();

  for (std::size_t i = 0; i < n_roots; i++) {
    TA::DistArray<Tile, Policy> guess(f_ia.world(), range, f_ia.shape());

    guess.fill(numeric_type(0.0));

    // fill in with 1.0
    //    std::size_t idx_i = i % n_a;
    //    std::size_t idx_a = i / n_a;

    std::size_t idx_i = n_i - 1;
    std::size_t idx_a = i;
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
