//
// Created by Chong Peng on 7/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EA_EOM_CCSD_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EA_EOM_CCSD_IMPL_H_

#include "mpqc/chemistry/qc/lcao/cc/eom/eom_preconditioner.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
void EOM_EA_CCSD<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  auto target_precision = ex_energy->target_precision(0);
  if (vector_threshold_ == 0.0) {
    vector_threshold_ = 10 * target_precision;
  }
  if (!this->computed()) {
    auto& world = this->wfn_world()->world();

    auto ccsd_energy =
        std::make_shared<Energy>(this->shared_from_this(), target_precision);
    // do CCSD energy
    CCSD<Tile, Policy>::evaluate(ccsd_energy.get());

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nEOM-EA-CCSD Electron Affinity\n";
    auto n_roots = ex_energy->n_roots();
    ExEnv::out0() << indent << "Number of roots: " << n_roots << "\n";
    ExEnv::out0() << indent << "Max number of vector per root: " << max_vector_
                  << "\n";
    ExEnv::out0() << indent
                  << "Threshold for norm of new vector: " << vector_threshold_
                  << "\n";

    // initialize guest vector
    ExEnv::out0() << indent << "\nInitialize Guess Vector \n";
    auto C = init_guess_vector(n_roots);

    // initialize intermediates
    ExEnv::out0() << indent << "\nInitialize Intermediates \n";
    auto imds = cc::compute_eom_intermediates(this->lcao_factory(),
                                              this->ao_factory(), this->t2(),
                                              this->t1(), this->is_df(), "ea");

    auto max_iter = this->max_iter_;
    auto result =
        ea_eom_ccsd_davidson_solver(C, imds, max_iter, target_precision);

    this->computed_ = true;
    ExcitationEnergy::Provider::set_value(
        ex_energy, std::vector<numeric_type>(result.data(),
                                             result.data() + result.size()));

    auto time1 = mpqc::fenced_now(world);
    ExEnv::out0() << "EOM-EA-CCSD Total Time: "
                  << mpqc::duration_in_s(time0, time1) << " S\n";
  }
}

template <typename Tile, typename Policy>
std::vector<typename EOM_EA_CCSD<Tile, Policy>::GuessVector>
EOM_EA_CCSD<Tile, Policy>::init_guess_vector(std::size_t n_roots) {
  auto& world = this->wfn_world()->world();

  std::vector<typename EOM_EA_CCSD<Tile, Policy>::GuessVector> guess_vector(
      n_roots);

  // occ range
  auto range_o = this->lcao_factory()
                     .orbital_registry()
                     .retrieve(OrbitalIndex(L"i"))
                     .trange();

  // unocc range
  auto range_v = this->lcao_factory()
                     .orbital_registry()
                     .retrieve(OrbitalIndex(L"a"))
                     .trange();

  TA::TiledRange trange_ca = {range_v};
  TA::TiledRange trange_cabi = {range_v, range_v, range_o};

#if TA_DEFAULT_POLICY == 1
  // making fake shape of 1
  TA::SparseShape<float> shape_v(TA::Tensor<float>(trange_ca.tiles_range(), 1),
                                 trange_ca);
  TA::SparseShape<float> shape_vvo(
      TA::Tensor<float>(trange_cabi.tiles_range(), 0.0), trange_cabi);
#endif

  for (std::size_t i = 0; i < n_roots; ++i) {
    // Ca
    guess_vector[i] = ::mpqc::cc::TPack<TArray>(2);

#if TA_DEFAULT_POLICY == 0
    guess_vector[i][0] = TA::DistArray<Tile, Policy>(world, trange_ca);
#elif TA_DEFAULT_POLICY == 1
    guess_vector[i][0] = TA::DistArray<Tile, Policy>(world, trange_ca, shape_v);
#endif

    guess_vector[i][0].fill(0.0);
    // make unit vector
    std::vector<std::size_t> element_idx = {i};

    auto tile_idx = trange_ca.element_to_tile(element_idx);

    if (guess_vector[i][0].is_local(tile_idx)) {
      auto tile = guess_vector[i][0].find(tile_idx).get();
      tile[element_idx] = numeric_type(1.0);
    }
    guess_vector[i][0].truncate();

    // Cabi
#if TA_DEFAULT_POLICY == 0
    guess_vector[i][1] = TA::DistArray<Tile, Policy>(world, trange_cabi);
#elif TA_DEFAULT_POLICY == 1
    guess_vector[i][1] =
        TA::DistArray<Tile, Policy>(world, trange_cabi, shape_vvo);
#endif
    guess_vector[i][1].fill(0.0);
    guess_vector[i][1].truncate();
  }
  return guess_vector;
}

template <typename Tile, typename Policy>
EigenVector<typename Tile::numeric_type>
EOM_EA_CCSD<Tile, Policy>::ea_eom_ccsd_davidson_solver(
    std::vector<typename EOM_EA_CCSD<Tile, Policy>::GuessVector>& C,
    const cc::Intermediates<TA::DistArray<Tile, Policy>>& imds,
    std::size_t max_iter, double convergence) {
  std::size_t n_roots = C.size();

  /// make preconditioner
  std::shared_ptr<DavidsonDiagPred<GuessVector>> pred;
  {
    EigenVector<numeric_type> eps_o =
        array_ops::array_to_eigen(imds.FIJ).diagonal();
    EigenVector<numeric_type> eps_v =
        array_ops::array_to_eigen(imds.FAB).diagonal();

    pred = std::make_shared<cc::EAPred<TArray>>(eps_o, eps_v);
  }

  /// make operator
  auto op = [this, &imds](const std::vector<GuessVector>& vec) {
    std::size_t dim = vec.size();
    //    ExEnv::out0() << "vector dimension: " << dim << std::endl;

    // compute product of H with guess vector
    std::vector<GuessVector> HC(dim);
    for (std::size_t i = 0; i < dim; ++i) {
      if (vec[i][0].is_initialized() && vec[i][1].is_initialized()) {
        HC[i] = ::mpqc::cc::TPack<TArray>(2);
        HC[i][0] = compute_HS1(vec[i][0], vec[i][1], imds);
        HC[i][1] = compute_HS2(vec[i][0], vec[i][1], imds);

      } else {
        throw ProgrammingError("Guess Vector not initialized", __FILE__,
                               __LINE__);
      }
    }

    return HC;

  };

  /// make davidson object
  DavidsonDiag<GuessVector> dvd(n_roots, false, 2, max_vector_,
                                vector_threshold_);

  EigenVector<numeric_type> eig = EigenVector<numeric_type>::Zero(n_roots);

  eig = dvd.solve(C, op, pred.get(), convergence, max_iter);

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, false);

  return eig;
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_EA_CCSD<Tile, Policy>::compute_HS1(
    const TArray& Ca, const TArray& Cabi,
    const cc::Intermediates<TA::DistArray<Tile, Policy>>& imds) {
  TArray HS1;

  {
    HS1("a") = imds.FAB("a,c") * Ca("c");

    HS1("a") +=
        imds.FIA("k,c") * (2.0 * Cabi("a,c,k") - Cabi("c,a,k")) +
        (2.0 * imds.Waibc("a,k,b,c") - imds.Waibc("a,k,c,b")) * Cabi("b,c,k");
  }
  return HS1;
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_EA_CCSD<Tile, Policy>::compute_HS2(
    const TArray& Ca, const TArray& Cabi,
    const cc::Intermediates<TA::DistArray<Tile, Policy>>& imds) {
  TArray HS2;
  {
    HS2("a,b,i") = imds.Wabci("a,b,c,i") * Ca("c");

    HS2("a,b,i") +=
        imds.FAB("a,c") * Cabi("c,b,i") + imds.FAB("b,c") * Cabi("a,c,i") -
        imds.FIJ("k,i") * Cabi("a,b,k") +
        (2 * imds.Wiabj("k,b,c,i") - imds.Wiajb("k,b,i,c")) * Cabi("a,c,k") -
        imds.Wiajb("k,a,i,c") * Cabi("c,b,k") -
        imds.Wiabj("k,b,c,i") * Cabi("c,a,k") +
        imds.Wabcd("a,b,c,d") * Cabi("c,d,i") -
        (2 * imds.Wijab("k,l,c,d") - imds.Wijab("k,l,d,c")) * Cabi("c,d,l") *
            this->t2()("a,b,k,i");
  }
  return HS2;
}

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EA_EOM_CCSD_IMPL_H_
