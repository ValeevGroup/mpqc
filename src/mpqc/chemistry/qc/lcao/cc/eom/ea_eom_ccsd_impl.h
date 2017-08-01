//
// Created by Chong Peng on 7/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EA_EOM_CCSD_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EA_EOM_CCSD_IMPL_H_

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
void EA_EOM_CCSD<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  auto target_precision = ex_energy->target_precision(0);
  if (!this->computed()) {
    auto& world = this->wfn_world()->world();

    auto ccsd_energy =
        std::make_shared<Energy>(this->shared_from_this(), target_precision);
    // do CCSD energy
    CCSD<Tile, Policy>::evaluate(ccsd_energy.get());

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nEOM-CCSD Electron Affinity\n";
    auto n_roots = ex_energy->n_roots();
    ExEnv::out0() << indent << "Number of roots: " << n_roots << "\n";
    ExEnv::out0() << indent << "Max number of vector per root: " << max_vector_
                  << "\n";
    ExEnv::out0() << indent
                  << "Threshold for norm of new vector: " << vector_threshold_
                  << "\n";

    // initialize guest vector
    ExEnv::out0() << indent << "\nInitialize Guess Vector \n";


    // initialize intermediates
    ExEnv::out0() << indent << "\nInitialize Intermediates \n";

    auto max_iter = this->max_iter_;
//    auto result = ea_eom_ccsd_davidson_solver(max_iter, target_precision);
    std::vector<double> result(n_roots);

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
std::vector<typename EA_EOM_CCSD<Tile, Policy>::GuessVector>
EA_EOM_CCSD<Tile, Policy>::init_guess_vector(std::size_t n_roots) {
  std::vector<typename EA_EOM_CCSD<Tile, Policy>::GuessVector> guess_vector(
      n_roots);

  for (std::size_t i = 0; i < n_roots; ++i) {
  }
}

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EA_EOM_CCSD_IMPL_H_
