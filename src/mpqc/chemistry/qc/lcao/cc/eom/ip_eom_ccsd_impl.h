//
// Created by Chong Peng on 7/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_IP_EOM_CCSD_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_IP_EOM_CCSD_IMPL_H_

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
void IP_EOM_CCSD<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  auto target_precision = ex_energy->target_precision(0);
  if (!this->computed()) {
    auto& world = this->wfn_world()->world();

    auto ccsd_energy =
        std::make_shared<Energy>(this->shared_from_this(), target_precision);
    // do CCSD energy
    CCSD<Tile, Policy>::evaluate(ccsd_energy.get());

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nIP-EOM-CCSD Ionization Potential\n";
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
    auto imds = compute_intermediates();

    auto max_iter = this->max_iter_;
    auto result =
        ip_eom_ccsd_davidson_solver(C, imds, max_iter, target_precision);

    this->computed_ = true;
    ExcitationEnergy::Provider::set_value(
        ex_energy, std::vector<numeric_type>(result.data(),
                                             result.data() + result.size()));

    auto time1 = mpqc::fenced_now(world);
    ExEnv::out0() << "IP-EOM-CCSD Total Time: "
                  << mpqc::duration_in_s(time0, time1) << " S\n";
  }
}

template <typename Tile, typename Policy>
std::vector<typename IP_EOM_CCSD<Tile, Policy>::GuessVector>
IP_EOM_CCSD<Tile, Policy>::init_guess_vector(std::size_t n_roots) {
  auto& world = this->wfn_world()->world();

  std::vector<typename IP_EOM_CCSD<Tile, Policy>::GuessVector> guess_vector(
      n_roots);

  // occ range
  auto range_o = this->lcao_factory()
                     .orbital_registry()
                     .retrieve(OrbitalIndex(L"i"))
                     .trange();
  std::size_t n_o = range_o.elements_range().second;

  // unocc range
  auto range_v = this->lcao_factory()
                     .orbital_registry()
                     .retrieve(OrbitalIndex(L"a"))
                     .trange();

  TA::TiledRange trange_ci = {range_o};
  TA::TiledRange trange_caij = {range_v, range_o, range_o};

  // making fake shape of 1
  TA::SparseShape<float> shape_o(TA::Tensor<float>(trange_ci.tiles_range(), 1),
                                 trange_ci);

  TA::SparseShape<float> shape_voo(
      TA::Tensor<float>(trange_caij.tiles_range(), 0.0), trange_caij);

  for (std::size_t i = 0; i < n_roots; ++i) {
    // Ci
    guess_vector[i].t1 = TA::DistArray<Tile, Policy>(world, trange_ci, shape_o);
    guess_vector[i].t1.fill(0.0);
    // make unit vector
    std::vector<std::size_t> element_idx{{n_o - i - 1}};

    auto tile_idx = trange_ci.element_to_tile(element_idx);

    if (guess_vector[i].t1.is_local(tile_idx)) {
      auto tile = guess_vector[i].t1.find(tile_idx).get();
      tile[element_idx] = numeric_type(1.0);
    }
    guess_vector[i].t1.truncate();

    // Caij
    guess_vector[i].t2 =
        TA::DistArray<Tile, Policy>(world, trange_caij, shape_voo);
    guess_vector[i].t2.fill(0.0);
    guess_vector[i].t2.truncate();
  }
  return guess_vector;
}

template <typename Tile, typename Policy>
cc::Intermediates<Tile, Policy>
IP_EOM_CCSD<Tile, Policy>::compute_intermediates() {
  bool df = this->is_df();

  return cc::compute_intermediates(this->lcao_factory(), this->ao_factory(),
                                   this->t2(), this->t1(), df, "ip");
}

template <typename Tile, typename Policy>
EigenVector<typename Tile::numeric_type>
IP_EOM_CCSD<Tile, Policy>::ip_eom_ccsd_davidson_solver(
    std::vector<typename IP_EOM_CCSD<Tile, Policy>::GuessVector>& C,
    const cc::Intermediates<Tile, Policy>& imds, std::size_t max_iter,
    double convergence) {
  madness::World& world =
      C[0].t1.is_initialized() ? C[0].t1.world() : C[0].t2.world();
  std::size_t iter = 0;
  std::size_t n_roots = C.size();
  double norm_r = 1.0;

  /// make preconditioner
  Preconditioner pred;
  {
    std::cout << imds.FIJ << std::endl;
    EigenVector<numeric_type> eps_o =
        array_ops::array_to_eigen(imds.FIJ).diagonal();
    EigenVector<numeric_type> eps_v =
        array_ops::array_to_eigen(imds.FAB).diagonal();

    pred = Preconditioner(eps_o, eps_v);
  }

  /// make davidson object
  DavidsonDiag<GuessVector> dvd(n_roots, false, 2, max_vector_,
                                vector_threshold_);

  EigenVector<numeric_type> eig = EigenVector<numeric_type>::Zero(n_roots);

  while (iter < max_iter && norm_r > convergence) {
    auto time0 = mpqc::fenced_now(world);
    std::size_t dim = C.size();
    //    ExEnv::out0() << "vector dimension: " << dim << std::endl;

    // compute product of H with guess vector
    std::vector<GuessVector> HC(dim);
    for (std::size_t i = 0; i < dim; ++i) {
      if (C[i].t1.is_initialized() && C[i].t2.is_initialized()) {
        HC[i].t1 = compute_HS1(C[i].t1, C[i].t2, imds);
        HC[i].t2 = compute_HS2(C[i].t1, C[i].t2, imds);

      } else {
        throw ProgrammingError("Guess Vector not initialized", __FILE__,
                               __LINE__);
      }
    }

    auto time1 = mpqc::fenced_now(world);
    EigenVector<double> eig_new = dvd.extrapolate(HC, C, pred);
    auto time2 = mpqc::fenced_now(world);

    EigenVector<numeric_type> delta_e = eig - eig_new;
    norm_r = delta_e.norm();

    util::print_excitation_energy_iteration(iter, delta_e, eig_new,
                                            mpqc::duration_in_s(time0, time1),
                                            mpqc::duration_in_s(time1, time2));

    eig = eig_new;
    iter++;

  }  // end of while loop

  if (iter == max_iter) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter, "IP-EOM-CCSD");
  }

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, false);

  return eig;
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> IP_EOM_CCSD<Tile, Policy>::compute_HS1(
    const TArray& Ci, const TArray& Caij,
    const cc::Intermediates<Tile, Policy>& imds) {
  TArray HS1;

  {
    HS1("i") = -imds.FIJ("k,i") * Ci("k");

    HS1("i") +=
        2 * imds.FIA("l,d") * Caij("d,i,l") -

        imds.FIA("k,d") * Caij("d,k,i") -

        (2 * imds.Wijka("k,l,i,d") - imds.Wijka("l,k,i,d")) * Caij("d,k,l");
  }

  return HS1;
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> IP_EOM_CCSD<Tile, Policy>::compute_HS2(
    const TArray& Ci, const TArray& Caij,
    const cc::Intermediates<Tile, Policy>& imds) {
  TArray HS2;
  {
    HS2("a,i,j") = -imds.Wiajk("k,a,i,j") * Ci("k");

    HS2("a,i,j") += imds.FAB("a,c") * Caij("c,i,j") -

                    imds.FIJ("k,i") * Caij("a,k,j") -

                    imds.FIJ("l,j") * Caij("a,i,l") +

                    imds.Wijkl("k,l,i,j") * Caij("a,k,l") +

                    2 * imds.Wiabj("l,a,c,j") * Caij("c,i,l") -

                    imds.Wiabj("k,a,c,j") * Caij("c,k,i") -

                    imds.Wiajb("l,a,j,c") * Caij("c,i,l") -

                    imds.Wiajb("k,a,i,c") * Caij("c,k,j") -

                    (2 * imds.Wijab("l,k,d,c") - imds.Wijab("k,l,d,c")) *
                        Caij("d,k,l") * this->t2()("c,a,i,j");
  }
  return HS2;
}

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_IP_EOM_CCSD_IMPL_H_
