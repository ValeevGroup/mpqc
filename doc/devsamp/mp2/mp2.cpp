#include "mpqc/chemistry/qc/lcao/scf/rhf.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

using namespace mpqc;

using LCAOWfn = lcao::LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;

// This is a basic implementation of (iterative) MP2 energy.
// This is formulated in terms of occupied and unoccupied states represented
// as LCAOs, hence it's derived from LCAOWfn, aka ::mpqc::lcao::LCAOWavefunction
// .
// This class can only compute the energy, this is indicated by deriving it from
// Provides<Energy> (this introduces virtual methods can_evaluate() and
// evaluate() that specify, respectively, the ability to compute the energy and
// how the energy is computed).
class MP2 : public LCAOWfn, public Provides<Energy> {
 public:
  // a few abbreviations to make the code less verbose
  using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;
  using RHF = lcao::RHF<TA::TensorD, TA::SparsePolicy>;

  // the KeyVal constructor takes a KeyVal object that represents a keyword
  // group of the input file that corresponds to this object (see the "mp2"
  // group in the mp2.json file that accompanies this example).
  // The Keyval object will be queried for all keywords needed by
  // the KeyVal ctor of LCAOWfn, as well as keyword "ref" that specifies
  // the reference wave function.
  MP2(const KeyVal& kv) : LCAOWfn(kv) {
    ref_wfn_ = kv.class_ptr<RHF>("ref");
    if (!ref_wfn_)
      throw InputError("missing reference RHF wave function", __FILE__,
                       __LINE__, "ref");
  }

 private:
  // This implements the Energy::Provider::can_evaluate() virtual
  // function. This function returns true is the energy object can be computed.
  bool can_evaluate(Energy* energy) override {
    // can only compute energies (not forces (energy->order() == 1),
    // hessians, or higher derivatives)
    return energy->order() == 0;
  }

  // This implements the Energy::Provider::evaluate() virtual function.
  // This function computes the MP2 energy and assigns it to the Energy object.
  void evaluate(Energy* energy) override {
    // how precisely to compute the energy
    auto target_precision = energy->target_precision(0);

    // if has not been computed to the desired precision
    // (or never computed at all) ...
    if (computed_precision_ > target_precision) {
      // compute reference to higher precision than this wfn
      auto target_ref_precision = target_precision / 100.;
      auto ref_energy =
          std::make_shared<Energy>(ref_wfn_, target_ref_precision);
      ref_wfn_->evaluate(ref_energy.get());

      // use the reference orbitals to populate the orbital space registry
      init_sdref(ref_wfn_, target_ref_precision);

      // this actually computes the energy
      double mp2_corr_energy = compute_mp2_energy(target_precision);

      energy_ = ref_energy->energy() + mp2_corr_energy;
    }

    // commit the result to energy
    this->set_value(energy, energy_);
  }

  // This function actually solves the MP1 equations and returns the MP2 energy.
  double compute_mp2_energy(double target_precision) {
    // fac is an LCAOFactory object which evaluates integrals in terms of AOs
    // and LCAOs
    auto& fac = this->lcao_factory();
    auto& world = fac.world();
    // ofac is an OrbitalSpaceRegistry that defines the LCAO spaces that fac can
    // use ofac was populated by the init_sdref() call above
    auto& ofac = fac.orbital_registry();

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

    // lambda function will be used to do a Jacobi update of the residual
    auto jacobi_update = [eps_o, eps_v](TA::TensorD& result_tile) {

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

    // solve the MP1 equations
    auto converged = false;
    auto iter = 0;
    auto energy = +1.0;
    ExEnv::out0() << "Start solving MP2 Energy\n";
    while (not converged) {
      Array R;
      R("i,a,j,b") = G("i,a,j,b") + Fv("a,c") * T_("i,c,j,b") +
                     Fv("b,c") * T_("i,a,j,c") - Fo("i,k") * T_("k,a,j,b") -
                     Fo("j,k") * T_("i,a,k,b");

      // estimate the MP2 energy ... the Hylleraas formula is quadratic in error
      double updated_energy =
          (G("i,a,j,b") + R("i,a,j,b")).dot(2 * T_("i,a,j,b") - T_("i,b,j,a"));
      ExEnv::out0() << indent << "Iteration: " << iter
                    << " Energy: " << updated_energy << std::endl;

      computed_precision_ = std::abs(updated_energy - energy);
      energy = updated_energy;

      // update the amplitudes, if needed
      converged = computed_precision_ <= target_precision;
      if (not converged) {
        // R^{ij}_{ab} -> R^{ij}_{ab} / (F^i_i + F^j_j - F^a_a - F^b_b)
        TA::foreach_inplace(R, jacobi_update);
        // need a fence here since foreach_inplace mutates the contents of R
        // as a side effect.
        // N.B. most TiledArray ops will not need a fence (except to control
        //      the resource use)
        world.gop.fence();
        T_("i,a,j,b") += R("i,a,j,b");
        ++iter;
      }
    }

    return energy;
  }

  // This reimplements the ::mpqc::Wavefunction::obsolete() virtual function.
  // It gets called when, for example, the atomic positions get updated in
  // geometry optimization.
  void obsolete() override {
    LCAOWfn::obsolete();
    ref_wfn_->obsolete();
    computed_precision_ = std::numeric_limits<double>::max();
    energy_ = 0.0;
  }

  std::shared_ptr<RHF> ref_wfn_;
  Array T_;
  double computed_precision_ = std::numeric_limits<double>::max();
  double energy_ = 0.0;
};

// This macro registers the KeyVal constructor of our MP2 class and associates
// it with the "MP2" key, so that the KeyVal class knows how to find it when it
// finds this key as the object type in the user input.
MPQC_CLASS_EXPORT2("MP2", MP2);

// Creating this variable forces the code for the MP2 class to be linked into
// the mp2 executable (otherwise the MPQC main function will not see any
// references to this class and thus the linker will simply skip it).
mpqc::detail::ForceLink<MP2> fl;
