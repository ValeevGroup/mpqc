//
//
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_CCSD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_CC_EOM_EOM_CCSD_H_

#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/chemistry/qc/lcao/cc/ccsd_hbar.h"
#include "mpqc/chemistry/qc/lcao/ci/cis.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/math/linalg/davidson_diag.h"

namespace mpqc {
namespace lcao {

// close-shell eom-ccsd
// still working on it
template <typename Tile, typename Policy>
class EOM_CCSD : public CCSD<Tile, Policy>, public Provides<ExcitationEnergy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using GuessVector = ::mpqc::cc::T1T2<TArray, TArray>;
  using numeric_type = typename Tile::numeric_type;

 private:

  // preconditioner in DavidsonDiag, approximate the diaagonal H_bar matrix
  struct Preconditioner {
    /// diagonal of F_ij matrix
    EigenVector<numeric_type> eps_o;
    /// diagonal of F_ab matrix
    EigenVector<numeric_type> eps_v;

    Preconditioner(const EigenVector<numeric_type> &eps_O,
                   const EigenVector<numeric_type> &eps_V)
        : eps_o(eps_O), eps_v(eps_V) {}

    // default constructor
    Preconditioner() : eps_o(), eps_v() {}

    void operator()(const numeric_type &e, GuessVector &guess) const {
      const auto &eps_v = this->eps_v;
      const auto &eps_o = this->eps_o;

      auto task1 = [&eps_v, &eps_o, e](Tile &result_tile) {
        const auto &range = result_tile.range();
        float norm = 0.0;
        for (const auto &i : range) {
          const auto result = result_tile[i] / (e + eps_o[i[1]] - eps_v[i[0]]);
          result_tile[i] = result;
          norm += result * result;
        }
        return std::sqrt(norm);
      };

      auto task2 = [&eps_v, &eps_o, e](Tile &result_tile) {
        const auto &range = result_tile.range();
        float norm = 0.0;
        for (const auto &i : range) {
          const auto result = result_tile[i] / (e - eps_v[i[0]] - eps_v[i[1]] +
                                                eps_o[i[2]] + eps_o[i[3]]);
          result_tile[i] = result;
          norm += result * result;
        }
        return std::sqrt(norm);
      };

      TA::foreach_inplace(guess.t1, task1);
      TA::foreach_inplace(guess.t2, task2);

      guess.t1.world().gop.fence();
    }
  };

  TArray Fab_;
  TArray Fij_;
  TArray Fai_;
  TArray Gijkl_;
  TArray Gijka_;
  TArray Gabij_;
  TArray Giajb_;
  TArray Giabc_;
  TArray Gaibc_;
  TArray Gabcd_;

  TArray F_;
  TArray FAB_;
  TArray FIJ_;
  TArray FIA_;

  TArray WIbAj_;
  TArray WIbaJ_;

  TArray WAbCd_;
  TArray WAbCi_;

  TArray WKlIj_;
  TArray WKaIj_;

  TArray WAkCd_;
  TArray WKlIc_;
  // TArray WKliC_;

  std::vector<GuessVector> C_;

  // compute F and W intermediates
  void compute_FWintermediates();

  // compute contractions of HSS, HSD, HDS, and HDD
  //                         with guess vector Ci
  // reference: CPL, 248 (1996), 189
  TArray compute_HSSC(TArray Cai);
  TArray compute_HSDC(TArray Cabij);
  TArray compute_HDSC(TArray Cai);
  TArray compute_HDDC(TArray Cabij);

  void init() {
    Fij_ = this->get_fock_ij();
    Fab_ = this->get_fock_ab();
    Fai_ = this->get_fock_ai();
    Gijkl_ = this->get_ijkl();
    Gijka_ = this->get_ijka();
    Gabij_ = this->get_abij();
    Giajb_ = this->get_iajb();
    Gaibc_ = this->get_aibc();
    Giabc_ = this->get_iabc();
    Gabcd_ = this->get_abcd();
    compute_FWintermediates();
  }

 public:
  EOM_CCSD(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {}
  // read guess vectors from input
  //    void read_guess_vectors(rapidjson::Document& in);

  void obsolete() override {
    CCSD<Tile, Policy>::obsolete();
    TArray F_ = TArray();
    TArray Gijkl_ = TArray();
    TArray Gijka_ = TArray();
    TArray Gabij_ = TArray();
    TArray Giajb_ = TArray();
    TArray Giabc_ = TArray();
    TArray Gaibc_ = TArray();
    TArray Gabcd_ = TArray();

    TArray FAB_ = TArray();
    TArray FIJ_ = TArray();
    TArray FIA_ = TArray();

    TArray WIbAj_ = TArray();
    TArray WIbaJ_ = TArray();

    TArray WAbCd_ = TArray();
    TArray WAbCi_ = TArray();

    TArray WKlIj_ = TArray();
    TArray WKaIj_ = TArray();

    TArray WAkCd_ = TArray();
    TArray WKlIc_ = TArray();
  }

 protected:
  bool can_evaluate(ExcitationEnergy *ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy *ex_energy) override;

 private:
  // not complete
  void davidson_solver(std::size_t max_iter, double convergence);

  // compute energies (not complete, now just test intermediates)
  double compute_energy(std::size_t max_iter, double convergence);
};

#if TA_DEFAULT_POLICY == 0
extern template class EOM_CCSD<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "eom_ccsd_impl.h"

#endif
