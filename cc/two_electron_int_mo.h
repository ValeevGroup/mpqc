//
// Created by Chong Peng on 7/6/15.
//

#ifndef TILECLUSTERCHEM_TWO_ELECTRON_INT_MO_H
#define TILECLUSTERCHEM_TWO_ELECTRON_INT_MO_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

namespace tcc{

  namespace cc {
    // class to compute and store the two electron integrals needed for ccsd

    template<typename Tile, typename Policy>
    class TwoElectronIntMO {
    public:

      typedef TA::Array <double, 2, Tile, Policy> TArray2;
      typedef TA::Array <double, 3, Tile, Policy> TArray3;
      typedef TA::Array <double, 4, Tile, Policy> TArray4;

      TwoElectronIntMO(const TArray3 &Xpq, const TArray2 &Ci,
                       const TArray2 &Ca) {
        {
          TArray3 Xab, Xij, Xai;
          Xab("X,a,b") = Xpq("X,mu,nu") * Ca("nu,a") * Ca("mu,b");
          Xij("X,i,j") = Xpq("X,mu,nu") * Ci("nu,i") * Ci("mu,j");
          Xai("X,a,i") = Xpq("X,mu,nu") * Ca("nu,a") * Ci("mu,i");
          // integrals

          // <ab|ij>
          abij_("a,b,i,j") = Xai("X,a,i") * Xai("X,b,j");
//        std::cout << abij_ << std::endl;
          // <ij|kl>
          ijkl_("i,j,k,l") = Xij("X,i,k") * Xij("X,j,l");
//        std::cout << ijkl_ << std::endl;
          // <ab|cd>
          abcd_("a,b,c,d") = Xab("X,a,c") * Xab("X,b,d");
          //std::cout << abcd_ << std::endl;
          // <ab|ci>
          abci_("a,b,c,i") = Xab("X,a,c") * Xai("X,b,i");
//        std::cout << abci_ << std::endl;
          // <ab|ic>
          abic_("a,b,i,c") = Xai("X,a,i") * Xab("X,b,c");
          // <ai|kl>
          aikl_("a,i,k,l") = Xai("X,a,k") * Xij("X,i,l");
//        std::cout << aikl_ << std::endl;
          // <ia|kl>
          iakl_("i,a,k,l") = Xai("X,a,l") * Xij("X,i,k");
          // <ai|bj>
          aibj_("a,i,b,j") = Xab("X,a,b") * Xij("X,i,j");
//        std::cout << aibj_ << std::endl;

        }
        TArray3::wait_for_lazy_cleanup(Xpq.get_world());
      }

      const TArray4 &get_abij() const {
        return abij_;
      }

      const TArray4 &get_ijkl() const {
        return ijkl_;
      }

      const TArray4 &get_abcd() const {
        return abcd_;
      }

      const TArray4 &get_abci() const {
        return abci_;
      }

      const TArray4 &get_abic() const {
        return abic_;
      }

      const TArray4 &get_aikl() const {
        return aikl_;
      }

      const TArray4 &get_iakl() const {
        return iakl_;
      }

      const TArray4 &get_aibj() const {
        return aibj_;
      }

      const TArray4 &get_aijb() const {
        return aijb_;
      }

    private:

      // two electron integrals
      // using notation <ab|ij>
      TArray4 abij_;
      TArray4 ijkl_;
      TArray4 abcd_;
      TArray4 abci_;
      TArray4 abic_;
      TArray4 aikl_;
      TArray4 iakl_;
      TArray4 aibj_;
      TArray4 aijb_;
    };

  }
}
#endif //TILECLUSTERCHEM_TWO_ELECTRON_INT_MO_H
