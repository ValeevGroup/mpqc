//
// Created by Chong Peng on 1/3/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_F12_DBGF2F12_H_
#define SRC_MPQC_CHEMISTRY_QC_F12_DBGF2F12_H_

#include "gf2f12.h"
#include "db_f12_intermediates.h"

namespace mpqc {
namespace lcao {

/**
 *  \brief Dual Basis GF2F12 class
 *  keyval name for this class DBGF2F12
 *
 */

template <typename Tile>
class DBGF2F12 : public GF2F12<Tile> {
 public:
  using typename GF2F12<Tile>::Policy;
  using typename GF2F12<Tile>::TArray;
  using typename GF2F12<Tile>::LCAOFactoryType;
  using typename GF2F12<Tile>::Matrix;
  using typename GF2F12<Tile>::real_t;

  DBGF2F12() = default;
  virtual ~DBGF2F12() = default;

  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords: all keywords from GF2F12
   */
  // clang-format on

  DBGF2F12(const KeyVal& kv) : GF2F12<Tile>(kv) {}

  using GF2F12<Tile>::value;

 private:

  /// override GF2F12's function to initialize obs and cabs orbitals
  void init() override {
    // initialize obs
    auto mol = this->wfn_world()->atoms();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
        this->lcao_factory(), orbital_energy, *mol, this->is_frozen_core(),
        this->occ_block(), this->unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

    if (this->use_cabs()) {
      // initialize cabs
      closed_shell_dualbasis_cabs_mo_build_svd(this->lcao_factory(),
                                               this->trange1_engine(), "VBS",
                                               this->unocc_block());
    }
  }

  /// override GF2F12's function to compute target orbital
  void init_target_orbital_diagonal() override {
    auto nfzc = this->trange1_engine()->get_nfrozen();
    auto nocc = this->trange1_engine()->get_active_occ();
    auto orbital = this->orbital();
    auto& world = this->wfn_world()->world();

    TArray C_x_ta;
    TA::TiledRange1 tr_x{0, 1};
    auto& orbital_registry = this->lcao_factory().orbital_space();
    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;

    if (orbital < 0) {
      auto o_space = orbital_registry.retrieve(OrbitalIndex(L"i"));
      auto C_o = array_ops::array_to_eigen(o_space.coefs());
      auto C_x = C_o.block(0, nocc + orbital, C_o.rows(), 1);
      auto tr_obs = o_space.coefs().trange().data().front();
      C_x_ta = array_ops::eigen_to_array<Tile, TA::SparsePolicy>(world, C_x,
                                                                 tr_obs, tr_x);

      auto x_space =
          OrbitalSpaceTArray(OrbitalIndex(L"x"), OrbitalIndex(L"κ"), C_x_ta);
      orbital_registry.add(x_space);
    } else if (orbital > 0) {
      auto v_space = orbital_registry.retrieve(OrbitalIndex(L"a"));
      auto C_v = array_ops::array_to_eigen(v_space.coefs());
      auto C_x = C_v.block(0, orbital-1, C_v.rows(), 1);
      auto tr_obs = v_space.coefs().trange().data().front();
      C_x_ta = array_ops::eigen_to_array<Tile, TA::SparsePolicy>(world, C_x,
                                                                 tr_obs, tr_x);

      auto x_space =
          OrbitalSpaceTArray(OrbitalIndex(L"x"), OrbitalIndex(L"Α"), C_x_ta);
      orbital_registry.add(x_space);
    }
  }

  /// override GF2F12's function to compute V term
  std::tuple<TArray, TArray> compute_V() override {
    return f12::compute_V_ixjy_ixyj_df(this->lcao_factory(), this->use_cabs());
  }
};

#if TA_DEFAULT_POLICY == 1
extern template class DBGF2F12<TA::TensorD>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_F12_DBGF2F12_H_
