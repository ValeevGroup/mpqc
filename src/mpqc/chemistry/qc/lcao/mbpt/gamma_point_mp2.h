#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_

#include "mpqc/chemistry/qc/lcao/mbpt/mp2.h"
#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"
#include "mpqc/chemistry/qc/lcao/mbpt/denom.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class GammaPointMP2 : public RMP2<Tile, Policy> {
public:
    using LCAOFactoryType = PeriodicLCAOFactory<Tile, Policy>;

    GammaPointMP2() = default;

    ~GammaPointMP2() = default;

    GammaPointMP2(const KeyVal &kv) : RMP2<Tile, Policy>(kv) {

        phf_wfn_ = kv.keyval("ref").class_ptr<PeriodicAOWavefunction<TA::TensorZ, Policy>>();
        lcao_factory_ = lcao::detail::construct_periodic_lcao_factory<Tile,Policy>(kv);
        print_detail_ = kv.value<bool>("print_detail", false);
    }

    void obsolete() override {
      e_mp2_ = 0.0;
      phf_wfn_->obsolete();
    }

private:
    double e_ref_;
    double e_mp2_;
    RowMatrixXd C_;
    Eigen::VectorXd eps_;
    std::shared_ptr<PeriodicAOWavefunction<TA::TensorZ, Policy>> phf_wfn_;

    bool print_detail_;
    std::shared_ptr<LCAOFactoryType> lcao_factory_;

private:
    bool can_evaluate(Energy *energy) override {
        // can only evaluate the energy
        return energy->order() == 0;
    }

    void evaluate(Energy *result) override {
        if (!this->computed()) {
            /// cast ref_wfn to Energy::Provider
            auto ref_evaluator = std::dynamic_pointer_cast<typename Energy::Provider>(phf_wfn_);
            if (ref_evaluator == nullptr) {
              std::ostringstream oss;
              oss << "RefWavefunction in GammaPointMP2" << phf_wfn_->class_key()
                  << " cannot compute Energy" << std::endl;
              throw InputError(oss.str().c_str(), __FILE__, __LINE__);
            }

            ref_evaluator->evaluate(result);

            double ref_energy = this->get_value(result).derivs(0)[0];

            auto time0 = mpqc::now(this->wfn_world()->world(), true);
            // initialize
            init();

            e_mp2_ = compute_gamma_point_mp2();
            auto time1 = mpqc::now(this->wfn_world()->world(), true);

            auto duration = mpqc::duration_in_s(time0, time1);
            if (print_detail_) {
              ExEnv::out0() << "\n Total Gamma-Point MP2 time: " << duration
                            << " s\n";
            }

            ExEnv::out0() << "\nGamma-Point MP2 Energy = " << e_mp2_ << std::endl;

            this->computed_ = true;
            this->set_value(result, ref_energy + e_mp2_);
        }
    }

    /*!
     * \brief This initializes gamma-point MP2
     */
    void init() override {
        auto nk = phf_wfn_->nk();
        auto k_size = 1 + detail::k_ord_idx(nk(0) - 1, nk(1) - 1, nk(2) - 1, nk);
        auto unitcell = lcao_factory_->pao_factory().unitcell();

        int64_t gamma_point;
        if (k_size % 2 == 1)
          gamma_point = (k_size - 1) / 2;
        else
          throw std::invalid_argument(
              "# of k points must be odd in order to run gamma-point methods");

        C_ = phf_wfn_->co_coeff()[gamma_point].real();
        eps_ = phf_wfn_->co_energy()[gamma_point].real();

        this->trange1_engine_ = mo_insert_gamma_point(*lcao_factory_, C_, unitcell, this->occ_block(),
                              this->unocc_block());

    }

    /*!
     * \brief This computes gamma-point MP2 energy
     * \return gamma-point MP2 energy
     */
    double compute_gamma_point_mp2() {
      auto &world = this->wfn_world()->world();

      ExEnv::out0() << "Computing conventional gamma-point MP2 ..." << std::endl;

      auto g_abij = lcao_factory_->compute(L"<a b|G|i j>");

      const auto charge = 0;
      const auto nelectrons =
          lcao_factory_->pao_factory().unitcell().total_atomic_number() -
          charge;
      auto n_occ = nelectrons / 2;
      std::size_t n_frozen = 0;
      auto t2 = d_abij(g_abij, eps_, n_occ, n_frozen);

      auto time0 = mpqc::now(world, false);

      double e_mp2 = TA::dot(
          2.0 * g_abij("a, b, i, j") - g_abij("b, a, i, j"), t2("a, b, i, j"));

      auto time1 = mpqc::now(world, false);
      auto duration = mpqc::duration_in_s(time0, time1);

      if (print_detail_) {
        ExEnv::out0() << " Energy computing time: " << duration << " s\n";
      }
      return e_mp2;
    }

};

#if TA_DEFAULT_POLICY == 0
extern template class GammaPointMP2<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class GammaPointMP2<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_GAMMA_POINT_MP2_H_
