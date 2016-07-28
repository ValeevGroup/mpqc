//
// Created by Chong Peng on 7/26/16.
//

#ifndef MPQC_CABS_SINGLES_H
#define MPQC_CABS_SINGLES_H

#include "../../../../../include/eigen.h"
#include "../../../../../include/tiledarray.h"
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>
#include <TiledArray/algebra/conjgrad.h>

namespace mpqc{
namespace f12{

template <typename Tile>
class CABSSingles {
public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;


  CABSSingles() = default;

  /**
   * @param mo_int MolecularIntegral Object
   * @param vir   if include F_ia in singles, default is true
   */
  CABSSingles(MolecularIntegralClass& mo_int, bool vir=true) : mo_int_(mo_int), couple_virtual_(vir){}

  real_t compute();

private:

  struct CABSSingleEquation {

    const TArray& F_AB_;
    const TArray& F_MN_;

    /**
     * @param F_AB  Fock matrix in all virtual space
     * @param F_MN  Fock matrix in occupied space
     */
    CABSSingleEquation(const TArray& F_AB, const TArray& F_MN)
            : F_AB_(F_AB), F_MN_(F_MN){}

    /**
     * @param[in] t T1 amplitude
     * @param[out] r  residual
     */
    void operator()(const TArray& t, TArray& r){
      r("i,A") = F_MN_("i,j")*t("j,A") - t("i,B")*F_AB_("B,A");
    }

  };

  /**
   *  compute the preconditioner X_i^A' from 1/(F_i^i - F_A'^A')
   * @param[out] P_MA  precoditioner 
   * @param[in] F_AB  Fock matrix in all virtual space
   * @param[in] F_MN  Fock matrix in occupied space
   */
  void compute_preconditioner(TArray &P_MA, const TArray &F_AB, const TArray &F_MN);

private:

  MolecularIntegralClass& mo_int_;
  bool couple_virtual_;
};

template <typename Tile>
typename CABSSingles<Tile>::real_t
CABSSingles<Tile>::compute() {

  bool df = mo_int_.atomic_integral().orbital_basis_registry()->have(OrbitalIndex(L"Κ"));

  TArray F_AB, F_MN;
  if(df){
    F_AB = mo_int_.compute(L"<A'|F|B'>[df]");
    F_MN = mo_int_.compute(L"<m|F|n>[df]");

  }
  else{
    F_AB = mo_int_.compute(L"<A'|F|B'>");
    F_MN = mo_int_.compute(L"<m|F|n>");
  }


  TArray F_MA;
  /// include contribution of F_m^a into F_m^A'
  if(couple_virtual_){
    if(df){
      F_MA = mo_int_.compute(L"<m|F|A'>[df]");
    }
    else{
      F_MA = mo_int_.compute(L"<m|F|A'>");
    }
  }
  /// not include contribution of F_m^a into F_m^A'
  else{
    TArray F_Ma;
    if(df){
      F_Ma = mo_int_.compute(L"<m|F|a'>[df]");
    }
    else{
      F_Ma = mo_int_.compute(L"<m|F|a'>");
    }
    MatrixD F_Ma_eigen = array_ops::array_to_eigen(F_Ma);
    auto n_occ = F_Ma_eigen.rows();
    auto n_cabs = F_Ma_eigen.cols();
    auto n_allvir = F_AB.trange().elements().extent()[0];
    auto n_vir = n_allvir - n_cabs;

    MatrixD F_MA_eigen = MatrixD::Zero(n_occ, n_allvir);
    F_MA_eigen.block(0,n_vir,n_occ,n_cabs) << F_Ma_eigen;

    auto tr_m = F_Ma.trange().data()[0];
    auto tr_A = F_AB.trange().data()[0];

    F_MA = array_ops::eigen_to_array<Tile>(F_Ma.get_world(), F_MA_eigen,tr_m, tr_A);
    F_MA.truncate();

  }
//  std::cout << F_MA << std::endl;

  TArray t;
  t("i,A") = -F_MA("i,A");

  // compute preconditioner
  TArray P_MA(F_MA.get_world(), F_MA.trange(), F_MA.get_shape());
  compute_preconditioner(P_MA, F_AB, F_MN);
//  std::cout << P_MA << std::endl;

  CABSSingleEquation cabs_singles(F_AB,F_MN);

  TA::ConjugateGradientSolver<TArray,CABSSingleEquation> cg_solver;
  auto resnorm = cg_solver(cabs_singles, F_MA, t, P_MA, 1e-12);


  real_t E_S = 2*TA::dot(t("i,A"), F_MA("i,A"));

  return E_S;

}

template <typename Tile>
void CABSSingles<Tile>::compute_preconditioner(TA::DistArray <Tile, TA::SparsePolicy> &P_MA,
                                          const TA::DistArray <Tile, TA::SparsePolicy> &F_AB,
                                          const TA::DistArray <Tile, TA::SparsePolicy> &F_MN)
{
  Eigen::MatrixXd F_AB_eigen = array_ops::array_to_eigen(F_AB);
  Eigen::MatrixXd F_MN_eigen = array_ops::array_to_eigen(F_MN);

  auto& world = P_MA.get_world();
  using range_type = typename TArray::range_type;
  using iterator = typename TArray::iterator;

  auto make_tile = [&F_AB_eigen, &F_MN_eigen](range_type &range){
    auto result_tile = Tile(range);
    // compute index
    const auto i0 = result_tile.range().lobound()[0];
    const auto in = result_tile.range().upbound()[0];
    const auto a0 = result_tile.range().lobound()[1];
    const auto an = result_tile.range().upbound()[1];

    auto ia = 0;
    typename Tile::value_type tmp = 1.0;
    for (auto i = i0; i < in; ++i) {
      const auto fi =  F_MN_eigen(i,i);
      for (auto a = a0; a < an; ++a, ++ia) {
        const auto fa = F_AB_eigen(a,a);
        const auto result_ai = tmp / (fi - fa);
        result_tile[ia] = result_ai;
      }
    }
    return result_tile;

  };

  for (iterator it = P_MA.begin(); it != P_MA.end(); ++it ){
    madness::Future<Tile> tile = world.taskq.add(
            make_tile,
            P_MA.trange().make_tile_range(it.ordinal())
            );
    *it = tile;
  }
  world.gop.fence();

}

} // end of namespace f12
} // end of namespace mpqc

#endif //MPQC_CABS_SINGLES_H
