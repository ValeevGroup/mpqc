#include "mpqc/chemistry/qc/lcao/scf/rhf.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

using namespace mpqc;

using LCAOWfn = lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;
using RHF = lcao::RHF<TA::TensorD, TA::SparsePolicy>;

/// basic (iterative) MP2 energy
class MP2 : public LCAOWfn, public Provides<Energy> {
 public:
  using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;

  /// the KeyVal ctor
  MP2(const KeyVal& kv) : LCAOWfn(kv) {
    ref_wfn_ = kv.class_ptr<RHF>("ref");
    if (!ref_wfn_)
      throw InputError("missing reference RHF wave function", __FILE__,
                       __LINE__, "ref");
  }

 private:
  /// this implements Energy::Provider::can_evaluate()
  bool can_evaluate(Energy* energy) override {
    // can only compute energies (not forces or hessians)
    return energy->order() == 0;
  }

  /// this implements Energy::Provider::evaluate()
  void evaluate(Energy* energy) override {

    // how precisely to compute the energy
    auto target_precision = energy->target_precision(0);

    // if has not been computed to the desired precision (or never computed at all) ...
    if (computed_precision_ > target_precision) {

      // compute reference to higher precision than this wfn
      auto target_ref_precision = target_precision / 100.;
      auto ref_energy =
          std::make_shared<Energy>(ref_wfn_, target_ref_precision);
      ref_wfn_->evaluate(ref_energy.get());

      /// use the reference orbitals to populate the orbital space registry
      init_sdref(ref_wfn_, target_ref_precision);

      /// compute
      double mp2_corr_energy = compute_mp2_energy(target_precision);

      /// commit the result
      this->set_value(energy, ref_energy->energy() + mp2_corr_energy);
      computed_precision_ = target_precision;
    }
  }

  /// this method actually solves the MP1 equations
  double compute_mp2_energy(double target_precision) {
    auto& fac = this->lcao_factory();
    auto& world = fac.world();
    auto& ofac = fac.orbital_space();

    auto nocc = ofac.retrieve("m").rank();
    auto nocc_act = ofac.retrieve("i").rank();
    auto nvir = ofac.retrieve("a").rank();
    auto nfzc = nocc - nocc_act;

    auto F = fac.compute(L"(p|F|q)");
    Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
    // replicated diagonal elements of Fo
    auto eps_o = eps_p.segment(nfzc, nocc_act);
    // replicated diagonal elements of Fv
    auto eps_v = eps_p.tail(nvir);

    // G_iajb
    auto G = fac.compute(L"(i a|G|j b)");
    // Fij
    auto Fo = fac.compute(L"(i|F|j)");
    // Fab
    auto Fv = fac.compute(L"(a|F|b)");

    // zero out amplitudes
    if (!T_.is_initialized()) {
      T_ = Array(world, G.trange(), G.shape());
      T_.fill(0.0);
    }

    // lambda function to update the residual
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
    auto converged = false;
    auto iter = 0;
    auto energy = +1.0;
    ExEnv::out0() << "Start solving MP2 Energy\n";
    while (not converged) {
      Array R;
      R("i,a,j,b") = G("i,a,j,b") + Fv("a,c") * T_("i,c,j,b") +
                     Fv("b,c") * T_("i,a,j,c") - Fo("i,k") * T_("k,a,j,b") -
                     Fo("j,k") * T_("i,a,k,b");

      double updated_energy =
          (G("i,a,j,b") + R("i,a,j,b")).dot(2 * T_("i,a,j,b") - T_("i,b,j,a"));
      ExEnv::out0() << indent << "Iteration: " << iter
                    << " Energy: " << updated_energy << std::endl;

      computed_precision_ = std::abs(updated_energy - energy);
      energy = updated_energy;

      // update the amplitudes, if needed
      converged = computed_precision_ <= target_precision;
      if (not converged) {
        TA::foreach_inplace(R, jacobi_update);
        world.gop
            .fence();  // <- need a fence for mutating ops like foreach_inplace
        T_("i,a,j,b") += R("i,a,j,b");
        ++iter;
      }
    }

    return energy;
  }

  void obsolete() override {
    computed_precision_ = std::numeric_limits<double>::max();
    LCAOWfn::obsolete();
  }

  std::shared_ptr<RHF> ref_wfn_;
  Array T_;
  double computed_precision_ = std::numeric_limits<double>::max();
};

MPQC_CLASS_EXPORT2("MP2", MP2);
mpqc::detail::ForceLink<MP2> fl;
