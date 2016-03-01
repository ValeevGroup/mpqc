//
//
#include "../include/tiledarray.h"
namespace mpqc{

  // close-shell eom-ccsd
  // still working on it
  class EOM_CCSD {
    //using CC = cc::CCSD<typename Tile, typename Policy>;
    //using CCinter = cc::CCSDIntermediate<typename Tile, typename Policy>;

    using TArray = TiledArray::TSpArrayD;

    TArray T1_;
    TArray T2_;
    
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
    TArray WKbIj_;

    TArray WAkCd_;
    TArray WKlIc_;
    TArray WKliC_;

    // compute F and W intermediates
    void compute_intermediates();

    public:
    template <typename CC, typename CCinter>
    EOM_CCSD(const CC& cc, CCinter& cc_inter, const TArray& fock) :
      T1_(cc.get_t1()),
      T2_(cc.get_t2()),
      F_(fock),
      Gijkl_(cc_inter->get_ijkl()),
      Gijka_(cc_inter->get_ijka()),
      Gabij_(cc_inter->get_abij()),
      Giajb_(cc_inter->get_iajb()),
      Giabc_(cc_inter->get_iabc()),
      Gaibc_(cc_inter->get_aibc()),
      Gabcd_(cc_inter->get_abcd())
    {
    }

    double compute_energy();

  };
}
