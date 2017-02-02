#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

using namespace mpqc;

class MP2 : public lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>,
            public CanEvaluate<Energy> {
 public:
  MP2(const KeyVal& kv)
      : lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.class_ptr<Wavefunction>("ref");
    } else {
      throw InputError("Default RefWavefunction in MP2 is not support! \n",
                       __FILE__, __LINE__, "ref");
    }
  }
  virtual ~MP2() = default;

 protected:
  bool can_evaluate(Energy* energy) override { return energy->order() == 0; }

  void evaluate(Energy* energy) override {
    if (!this->computed()) {
      /// cast ref_wfn to Energy::Evaluator
      auto ref_evaluator =
          std::dynamic_pointer_cast<typename Energy::Evaluator>(ref_wfn_);
      if (ref_evaluator == nullptr) {
        std::ostringstream oss;
        oss << "RefWavefunction in MP2" << ref_wfn_->class_key()
            << " cannot compute Energy" << std::endl;
        throw InputError(oss.str().c_str(), __FILE__, __LINE__, "ref");
      }

      ref_evaluator->evaluate(energy);

      double ref_energy = this->get_value(energy).derivs(0)[0];

      // initialize
      this->init();

      // compute
      compute_mp2(energy->target_precision(0));

      this->computed_ = true;
      this->set_value(energy, ref_energy + mp2_corr_energy_);
    }
  }

 private:
  void compute_mp2(double precision) {
    auto& world = this->wfn_world()->world();

    const auto n_active_occ = this->trange1_engine()->get_active_occ();
    const auto n_frozen = this->trange1_engine()->get_nfrozen();
    const auto n_vir = this->trange1_engine()->get_vir();

    const auto& ens = *this->orbital_energy();
    // replicated diagonal elements of Fo
    const auto eps_o = ens.segment(n_frozen, n_active_occ);
    // replicated diagonal elements of Fv
    const auto eps_v = ens.segment(n_frozen + n_active_occ, n_vir);

    // compute integrals
    auto& lcao_factory = this->lcao_factory();
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
    ExEnv::out0() << "Start solving MP2 Energy\n";
    while (not converged) {
      iter++;

      TA::DistArray<TA::TensorD, TA::SparsePolicy> R;

      R("i,a,j,b") = G("i,a,j,b") + Fv("a,c") * T("i,c,j,b") +
                     Fv("b,c") * T("i,a,j,c") - Fo("i,k") * T("k,a,j,b") -
                     Fo("j,k") * T("i,a,k,b");

      double norm = R("i,a,j,b").norm();
      ExEnv::out0() << "Iteration: " << iter << " Norm: " << norm << "\n";

      converged = norm < precision;

      if (not converged) {
        // update residual
        TA::foreach_inplace(R, jacobi_update);
        world.gop.fence();

        // update amplitudes
        T("i,a,j,b") += R("i,a,j,b");
      }
    }

    // compute energy
    mp2_corr_energy_ = (2 * G("i,a,j,b") - G("i,b,j,a")).dot(T("i,a,j,b"));
    ExEnv::out0() << "MP2 Correlation Energy = " << mp2_corr_energy_ << "\n";
  }

  // use default init() in LCAOWavefunction
  void init() override {
    LCAOWavefunction<TA::TensorD, TA::SparsePolicy>::init();
  }

 private:
  std::shared_ptr<lcao::Wavefunction> ref_wfn_;
  double mp2_corr_energy_;
};

MPQC_CLASS_EXPORT2("MP2", MP2);
mpqc::detail::ForceLink<MP2> fl;
