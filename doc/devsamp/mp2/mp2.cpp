#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/lcao/scf/rhf.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

using namespace mpqc;

using LCAOWfn = lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;
using RHF = lcao::RHF<TA::TensorD, TA::SparsePolicy>;

class MP2 : public LCAOWfn, public Provides<Energy> {
 public:
  MP2(const KeyVal& kv) : LCAOWfn(kv) {
    ref_wfn_ = kv.class_ptr<RHF>("ref");
    if (!ref_wfn_)
      throw InputError("missing reference RHF wave function", __FILE__, __LINE__,
                       "ref");
  }

 private:

  bool can_evaluate(Energy* energy) override { return energy->order() == 0; }

  void evaluate(Energy* energy) override {

    auto target_precision = energy->target_precision(0);
    // compute only if never computed, or requested with higher precision than before
    if (!this->computed() || computed_precision_ > target_precision) {

      // compute reference to higher precision than required of correlation energy
      auto target_ref_precision = target_precision / 100.;
      auto ref_energy = std::make_shared<Energy>(ref_wfn_, target_ref_precision);
      ref_wfn_->evaluate(ref_energy.get());

      init_sdref(ref_wfn_, target_ref_precision);

      // compute
      compute_mp2_energy(target_precision);

      this->computed_ = true;
      this->set_value(energy, ref_energy->energy() + mp2_corr_energy_);
      computed_precision_ = target_precision;
    }
  }

  void compute_mp2_energy(double precision) {
    auto& world = this->wfn_world()->world();

    auto& lcao_factory = this->lcao_factory();

    const auto n_active_occ = this->trange1_engine()->get_active_occ();
    const auto n_frozen = this->trange1_engine()->get_nfrozen();
    const auto n_vir = this->trange1_engine()->get_vir();

    auto F = lcao_factory.compute(L"(p|F|q)");
    Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
    // replicated diagonal elements of Fo
    const auto eps_o = eps_p.segment(n_frozen, n_active_occ);
    // replicated diagonal elements of Fv
    const auto eps_v = eps_p.segment(n_frozen + n_active_occ, n_vir);

    // compute integrals
    // G_iajb
    auto G = lcao_factory.compute(L"(i a|G|j b)");

    // Fij
    auto Fo = lcao_factory.compute(L"(i|F|j)");

    // Fab
    auto Fv = lcao_factory.compute(L"(a|F|b)");

    // T_iajb initialize as 0
    TA::DistArray<TA::TensorD, TA::SparsePolicy> T(world, G.trange(),
                                                   G.shape());
    T.fill(0.0);
    T.truncate();

    world.gop.fence();

    // lambda function to update residual
    auto jacobi_update = [&eps_o, &eps_v](TA::TensorD& result_tile) {

      const auto& range = result_tile.range();
      double norm = 0.0;
      for (const auto& i : range) {
        const auto result_abij = result_tile[i] / (eps_o[i[0]] - eps_v[i[1]] +
                                                   eps_o[i[2]] - eps_v[i[3]]);
        result_tile[i] = result_abij;
        norm += result_abij * result_abij;
      }
      return std::sqrt(norm);
    };

    // start iteration
    bool converged = false;
    std::size_t iter = 0;
    mp2_corr_energy_ = +1.0;
    ExEnv::out0() << "Start solving MP2 Energy\n";
    while (not converged) {
      iter++;

      TA::DistArray<TA::TensorD, TA::SparsePolicy> R;

      R("i,a,j,b") = G("i,a,j,b") + Fv("a,c") * T("i,c,j,b") +
                     Fv("b,c") * T("i,a,j,c") - Fo("i,k") * T("k,a,j,b") -
                     Fo("j,k") * T("i,a,k,b");

      double hylleraas_energy = (G("i,a,j,b") + R("i,a,j,b")).dot(2*T("i,a,j,b") - T("i,b,j,a"));
      ExEnv::out0() << indent << "Iteration: " << iter << " Energy: " << hylleraas_energy << "\n";

      const auto delta_energy = hylleraas_energy - mp2_corr_energy_;
      converged = std::abs(delta_energy) < precision;
      mp2_corr_energy_ = hylleraas_energy;

      if (not converged) {
        // update residual
        TA::foreach_inplace(R, jacobi_update);
        world.gop.fence();

        // update amplitudes
        T("i,a,j,b") += R("i,a,j,b");
      }
    }
  }

  std::shared_ptr<RHF> ref_wfn_;
  double mp2_corr_energy_;
  double computed_precision_ = std::numeric_limits<double>::max();
};

MPQC_CLASS_EXPORT2("MP2", MP2);
mpqc::detail::ForceLink<MP2> fl;
