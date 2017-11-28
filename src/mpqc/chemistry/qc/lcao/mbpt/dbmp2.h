//
// Created by Chong Peng on 6/13/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_H_

#include "mpqc/chemistry/qc/lcao/mbpt/mp2.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

namespace detail {

template <typename Tile>
struct ScfCorrection {
  using result_type = double;
  using argument_type = Tile;

  std::shared_ptr<Eigen::VectorXd> vec_;
  unsigned int n_occ_;

  ScfCorrection(std::shared_ptr<Eigen::VectorXd> vec, int n_occ)
      : vec_(std::move(vec)), n_occ_(n_occ) {}

  ScfCorrection(ScfCorrection const &) = default;

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

    auto sti = st[0];
    auto fni = fn[0];
    auto sta = st[1];
    auto fna = fn[1];

    for (auto i = sti; i < fni; ++i) {
      const auto e_i = vec[i];
      for (auto a = sta; a < fna; ++a, ++tile_idx) {
        const auto e_ia = e_i - vec[a + n_occ_];
        const auto data = tile.data()[tile_idx];
        me += (data * data) / (e_ia);
      }
    }
  }
};

template <typename Tile, typename Policy>
std::shared_ptr<::mpqc::utility::TRange1Engine> closed_shell_dual_basis_mo_build_steele(
    lcao::LCAOFactory<Tile, Policy> &lcao_factory,
    Eigen::VectorXd &ens,
    const Molecule &mols,
    bool frozen_core,
    std::size_t occ_blocksize,
    std::size_t vir_blocksize);

}  // namespace detail

/**
 *  \breif Dual basis MP2 method for closed shell system
 *
 * KeyVal type keyword: DBRMP2
 *
 */

template <typename Tile, typename Policy>
class DBRMP2 : public RMP2<Tile,Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = lcao::LCAOFactory<Tile, Policy>;

  DBRMP2() = default;

  virtual ~DBRMP2() { }

  DBRMP2(const KeyVal &kv) : RMP2<Tile,Policy>(kv) {
    method_ = kv.value<std::string>("method", "valeev");
    if( (method_!= "valeev") && (method_!="steele")){
      throw InputError("Invalid Method for Dual Basis MP2! \n",__FILE__,__LINE__,"method");
    }
  }

  void obsolete() override;

  /// function
  virtual double compute_scf_correction();

protected:

  // override RMP2's evaluate function
  void evaluate(Energy* result) override;

  // override RMP2's init
  void init();

private:
  std::string method_;
  double scf_correction_;
};

/**
 *  \breif Dual basis MP2 method for closed shell system with density-fitting
 *
 * KeyVal type keyword: RI-DBRMP2
 *
 */
template <typename Tile, typename Policy>
class RIDBRMP2 : public DBRMP2<Tile,Policy>{
public:
  /**
   * KeyVal constructor
   *
   * this class inherit all keywords from DBRMP2
   */

  RIDBRMP2(const KeyVal& kv) : DBRMP2<Tile,Policy>(kv) {};
  ~RIDBRMP2() { }

  /// override DBRMP2's compute_scf_correction function
  double compute_scf_correction() override;

private:

  /// override RMP2's compute function
  double compute() override;

};

#if TA_DEFAULT_POLICY == 0
extern template class DBRMP2<TA::TensorD, TA::DensePolicy>;
extern template class RIDBRMP2<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class DBRMP2<TA::TensorD, TA::SparsePolicy>;
extern template class RIDBRMP2<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "dbmp2_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_H_
