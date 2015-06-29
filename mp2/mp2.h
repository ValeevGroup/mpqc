//
// Created by Chong Peng on 6/24/15.
//

#ifndef TILECLUSTERCHEM_MP2_H
#define TILECLUSTERCHEM_MP2_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../ta_routines/array_to_eigen.h"
#include "trange1_engine.h"

using namespace tcc;

template <typename Tile, typename Policy>
class MP2{

public:

  typedef  TA::Array<double, 2, Tile, Policy> TArray2;
  typedef  TA::Array<double, 3, Tile, Policy> TArray3;
  typedef  TA::Array<double, 4, Tile, Policy> TArray4;


  MP2(){};

  //MP2(const TArray2 &fock, const TArray2 &s_ab, const TArray4 &abcd);

  MP2(const TArray2 &fock, const TArray2 &s_ab, const TArray3 &Xab,
        const TRange1Engine& tre) : tre_(tre)
  {
    auto tr1 = tre.get_all_tr1();
    init(fock, s_ab, Xab, tr1, tr1);
  };

  TArray4 get_two_e() const {return two_e_int_mo_;}

  void compute(){

    std::size_t occ_blocks = tre_.get_occ_blocks();
    std::size_t vir_blocks = tre_.get_vir_blocks();
    std::size_t occ_b = 0ul;
    std::size_t occ_e = occ_blocks;
    std::size_t vir_b = occ_blocks;
    std::size_t vir_e = occ_blocks + vir_blocks;
    std::vector<std::size_t> low {occ_b,vir_b,occ_b,vir_b};
    std::vector<std::size_t> up {occ_e,vir_e,occ_e,vir_e};
    TArray4 two_e_iajb;
    //two_e_iajb("i,a,j,b") = two_e_int_mo_("i,a,j,b").block({occ_b,vir_b,occ_b,vir_b}, {occ_e,vir_e,occ_e,vir_e});
    two_e_iajb("i,a,j,b") = two_e_int_mo_("i,a,j,b").block(low, up);
    double energy_mp2 = (two_e_iajb("i,a,j,b")*(2*two_e_iajb("i,a,j,b")-two_e_iajb("i,b,j,a")))
            .reduce(Mp2Red(en_mo_, tre_.get_occ()));

    //std::cout << two_e_iajb << std::endl;
    std::cout << energy_mp2 << std::endl;
  }

private:

  struct Mp2Red {
    using result_type = double;
    using argument_type = Tile;

    std::shared_ptr<Eig::VectorXd> vec_;
    unsigned int n_occ_;

    Mp2Red(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
            : vec_(std::move(vec)), n_occ_(n_occ) {}
    Mp2Red(Mp2Red const &) = default;

    result_type operator()() const { return 0.0; }
    result_type operator()(result_type const &t) const { return t; }
    void operator()(result_type &me, result_type const &other) const {
      me += other;
    }

    void operator()(result_type &me, argument_type const &tile) const {
      auto const &range = tile.range();
      auto const &vec = *vec_;
      auto const st = range.lobound_data();
      auto const fn = range.upbound_data();
      auto tile_idx = 0;
      for (auto i = st[0]; i < fn[0]; ++i) {
        const auto e_i = vec[i];
        for (auto a = st[1]; a < fn[1]; ++a) {
          const auto e_ia = e_i - vec[a + n_occ_];
          for (auto j = st[2]; j < fn[2]; ++j) {
            const auto e_iaj = e_ia + vec[j];
            for (auto b = st[3]; b < fn[3]; ++b, ++tile_idx) {
              const auto e_iajb = e_iaj - vec[b + n_occ_];
              me += 1 / (e_iajb)*tile.data()[tile_idx];
            }
          }
        }
      }
    }
  };

  //template <typename Tile, typename Policy>
  //void init(const TArray2& fock, const TArray2& s_mn, const TArray4& mnkl);

 // template <typename Tile, typename Policy>
  void init(const TArray2& fock, const TArray2& s_mn, const TArray3& Xmn,
            const TA::TiledRange1 &tr1, const TA::TiledRange1 &tr2 )
  {

    auto fock_eig = tcc::array_ops::array_to_eigen(fock);
    auto s_mn_eig = tcc::array_ops::array_to_eigen(s_mn);
    Eigen::GeneralizedSelfAdjointEigenSolver<decltype(s_mn_eig)> es(fock_eig, s_mn_eig);
    Eigen::VectorXd evals = es.eigenvalues();
    auto coeff_eig = es.eigenvectors();
    auto tr0 = Xmn.trange().data().back();
    coeff_mo_ = tcc::array_ops::eigen_to_array<Tile>(fock.get_world(),coeff_eig, tr0, tr1);

    //std::size_t col = tr0.elements().second;
    //std::size_t row = tr1.elements().second;
    //auto I = Eigen::MatrixXd::Identity(col, row);

    //auto I_TA = tcc::array_ops::eigen_to_array<Tile>(fock.get_world(), I, tr0, tr1);
    //I_TA.truncate();

    //reblocking Xmn
    //std::cout << Xmn << std::endl;
    //TArray3 Xmn_reblock;
    //Xmn_reblock("X,n,m1") = Xmn("X,m1,n1")*I_TA("m1,m");
    //std::cout << Xmn_reblock << std::endl;
    //Xmn_reblock("X,m,n") = Xmn_reblock("X,n,m1")*I_TA("m1,m");

    //std::cout << Xmn_reblock << std::endl;

    // construct two electron mo
    TArray3 Xmn_mo;
    Xmn_mo("X,m,n") = Xmn("X,mu,nu")*coeff_mo_("mu,m")*coeff_mo_("nu,n");
    two_e_int_mo_("p,q,r,s") = Xmn_mo("X,p,q")*Xmn_mo("X,r,s");
    en_mo_ = std::make_shared<Eigen::VectorXd>(std::move(evals));
  };


private:
  TArray4 two_e_int_mo_;
  TArray2 coeff_mo_;
  std::shared_ptr<Eigen::VectorXd> en_mo_;
  tcc::TRange1Engine tre_;
};

#endif //TILECLUSTERCHEM_MP2_H
