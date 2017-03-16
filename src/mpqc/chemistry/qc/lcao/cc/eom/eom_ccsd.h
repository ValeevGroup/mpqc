//
//
#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/chemistry/qc/lcao/ci/cis.h"
#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/math/linalg/davidson_diag.h"

namespace mpqc {
namespace lcao {

struct guess_vector {
  using TArray = TiledArray::TSpArrayD;
  TArray Cai;
  TArray Cabij;
};

// close-shell eom-ccsd
// still working on it
class EOM_CCSD : public CCSD<TA::TensorD, TA::SparsePolicy>,
                 public Provides<ExcitationEnergy> {
  // using CC = CCSD<typename Tile, typename Policy>;
  // using CCinter = CCSDIntermediate<typename Tile, typename Policy>;

  using TArray = TiledArray::TSpArrayD;

  TArray F_;
  TArray Gijkl_;
  TArray Gijka_;
  TArray Gabij_;
  TArray Giajb_;
  TArray Giabc_;
  TArray Gaibc_;
  TArray Gabcd_;

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

  std::vector<guess_vector> C_;

  // compute F and W intermediates
  void compute_FWintermediates();

  // compute contractions of HSS, HSD, HDS, and HDD
  //                         with guess vector Ci
  // reference: CPL, 248 (1996), 189
  TArray compute_HSSC(TArray Cai);
  TArray compute_HSDC(TArray Cabij);
  TArray compute_HDSC(TArray Cai);
  TArray compute_HDDC(TArray Cabij);

  void compute_HC();

  void init() {
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
  EOM_CCSD(const KeyVal &kv) : CCSD<TA::TensorD, TA::SparsePolicy>(kv) {}
  // read guess vectors from input
  //    void read_guess_vectors(rapidjson::Document& in);

  void obsolete() override {
    CCSD<TA::TensorD,TA::SparsePolicy>::obsolete();
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
    // TArray WKliC_ = TArray();
  }

  // not complete
  void davidson_solver(std::size_t max_iter, double convergence);

  // compute energies (not complete, now just test intermediates)
  double compute_energy(std::size_t max_iter, double convergence);

 protected:
  bool can_evaluate(ExcitationEnergy *ex_energy) override {
    return ex_energy->order() == 0;
  }

  void evaluate(ExcitationEnergy *ex_energy) override;
};

}  // namespace lcao
}  // namespace mpqc
