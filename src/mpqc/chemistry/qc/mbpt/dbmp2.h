//
// Created by Chong Peng on 6/13/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_H_

#include "mpqc/chemistry/qc/mbpt/mp2.h"

namespace mpqc {
namespace mbpt {

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
std::shared_ptr<TRange1Engine> closed_shell_dual_basis_mo_build_steele(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory,
    Eigen::VectorXd &ens,
    const Molecule &mols,
    bool frozen_core,
    std::size_t occ_blocksize,
    std::size_t vir_blocksize);

}  // end of namespace detail

/**
 * Dual basis MP2 method
 */

template <typename Tile, typename Policy>
class DBRMP2 : public RMP2<Tile,Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = integrals::LCAOFactory<Tile, Policy>;

  DBRMP2() = default;

  virtual ~DBRMP2();

  DBRMP2(const KeyVal &kv) : RMP2<Tile,Policy>(kv) {
    method_ = kv.value<std::string>("method", "valeev");
    if( (method_!= "valeev") && (method_!="steele")){
      throw std::runtime_error("Invalid Method for Dual Basis MP2! \n");
    }
  }

  double value() override;

  virtual double compute_scf_correction();

protected:
  void init() override;

private:
  std::string method_;

  //  std::shared_ptr<TRange1Engine> pulay_build(
  //      integrals::LCAOFactory<Tile, Policy> &lcao_factory,
  //      OrbitalSpaceRegistry<TArray>& orbital_registry,
  //      Eigen::VectorXd &ens, const rapidjson::Document &in,
  //      const Molecule &mols, int occ) {
  //    auto &ao_factory = lcao_factory.ao_factory();
  //    auto &world = ao_factory.world();
  //    using TArray = TA::DistArray<Tile, Policy>;
  //
  //    auto mo_time0 = mpqc::fenced_now(world);
  //    utility::print_par(world, "\nBuilding ClosedShell Dual Basis MO
  //    Orbital\n");
  //
  //    // solving occupied orbitals
  //    TArray F;
  //    if (ao_factory.registry().have(Formula(L"<μ|F|ν>"))) {
  //      F = ao_factory.registry().retrieve(Formula(L"<μ|F|ν>"));
  //    } else {
  //      F = ao_factory.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
  //    }
  //
  //    auto S = ao_factory.compute(L"<κ|λ>");
  //
  //    RowMatrixXd F_eig = array_ops::array_to_eigen(F);
  //    RowMatrixXd S_eig = array_ops::array_to_eigen(S);
  //
  //    // solve mo coefficients
  //    Eigen::GeneralizedSelfAdjointEigenSolver<RowMatrixXd> es(F_eig, S_eig);
  //
  //    bool frozen_core =
  //        in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
  //    std::size_t n_frozen_core = 0;
  //    if (frozen_core) {
  //      n_frozen_core = mols.core_electrons();
  //      utility::print_par(world, "Frozen Core: ", n_frozen_core, "
  //      electrons",
  //                         "\n");
  //      n_frozen_core = n_frozen_core / 2;
  //    }
  //
  //    Eigen::VectorXd ens_occ =
  //        es.eigenvalues().segment(n_frozen_core, occ - n_frozen_core);
  //    std::cout << es.eigenvalues() << std::endl;
  //    RowMatrixXd C_all = es.eigenvectors();
  //    RowMatrixXd C_occ = C_all.block(0, 0, S_eig.rows(), occ);
  //    RowMatrixXd C_corr_occ =
  //        C_all.block(0, n_frozen_core, S_eig.rows(), occ - n_frozen_core);
  //
  //    // finished solving occupied orbitals
  //
  //    // get all the sizes
  //    auto n_obs = S.trange().elements_range().extent()[0];
  //    std::size_t mo_blocksize =
  //        in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
  //    std::size_t vir_blocksize = in.HasMember("VirBlockSize")
  //                                    ? in["VirBlockSize"].GetInt()
  //                                    : mo_blocksize;
  //    std::size_t occ_blocksize = in.HasMember("OccBlockSize")
  //                                    ? in["OccBlockSize"].GetInt()
  //                                    : mo_blocksize;
  //    utility::print_par(world, "OccBlockSize: ", occ_blocksize, "\n");
  //    utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");
  //
  //    auto tre = std::make_shared<TRange1Engine>(occ, n_obs, occ_blocksize,
  //                                               vir_blocksize,
  //                                               n_frozen_core);
  //    auto tr_obs = S.trange().data().back();
  //    auto tr_occ = tre->compute_range(occ, occ_blocksize);
  //    auto tr_corr_occ = tre->get_active_occ_tr1();
  //    detail::parallel_print_range_info(world, tr_occ, "Occ");
  //    detail::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
  //
  //    // convert to TA
  //    auto C_occ_ta =
  //        array_ops::eigen_to_array<Tile>(world, C_occ, tr_obs, tr_occ);
  //    auto C_corr_occ_ta =
  //        array_ops::eigen_to_array<Tile>(world, C_corr_occ, tr_obs,
  //        tr_corr_occ);
  //
  //    // insert to registry
  //    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  //    auto occ_space =
  //        OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"),
  //        C_occ_ta);
  //    orbital_registry.add(occ_space);
  //
  //    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"),
  //                                             OrbitalIndex(L"κ"),
  //                                             C_corr_occ_ta);
  //    orbital_registry.add(corr_occ_space);
  //
  //    // solving the virtual orbitals
  //    // find fock matrix
  //    TArray F_vbs;
  //    // if use density fitting
  //    if (ao_factory.registry().have(Formula(L"<μ|F|ν>[df]"))) {
  //      F_vbs = ao_factory.compute(Formula(L"<Α|F|Β>[df]"));
  //    } else {
  //      F_vbs = ao_factory.compute(Formula(L"<Α|F|Β>"));
  //    }
  //    auto n_vbs = F_vbs.trange().elements_range().extent()[0];
  //    tre = std::make_shared<TRange1Engine>(occ, n_vbs, occ_blocksize,
  //                                          vir_blocksize, n_frozen_core);
  //
  //    // project Fock matrix against occ orbitals
  //    RowMatrixXd C_vir;
  //    {
  //      auto tr_vbs = F_vbs.trange().data().back();
  //
  //      // get identity matrix
  //      Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(n_vbs, n_vbs);
  //      auto identity_ta =
  //          array_ops::eigen_to_array<Tile>(world, identity, tr_vbs, tr_vbs);
  //      // get overlap
  //      auto S_vbs = ao_factory.compute(L"<Α|Β>");
  //      // get density
  //      RowMatrixXd D_obs = C_occ * C_occ.transpose();
  //
  //      RowMatrixXd D_vbs = RowMatrixXd::Zero(n_vbs, n_vbs);
  //      D_vbs.block(0, 0, n_obs, n_obs) << D_obs;
  //      TArray D_vbs_ta =
  //          array_ops::eigen_to_array<Tile>(world, D_vbs, tr_vbs, tr_vbs);
  //
  //      F_vbs("A,B") =
  //          (identity_ta("A,D") - 0.5 * S_vbs("A,mu") * D_vbs_ta("mu,D")) *
  //          F_vbs("D,E") *
  //          (identity_ta("E,B") - 0.5 * D_vbs_ta("E,mu") * S_vbs("mu,B"));
  //
  //      RowMatrixXd F_vbs_eigen = array_ops::array_to_eigen(F_vbs);
  //
  //      Eigen::SelfAdjointEigenSolver<RowMatrixXd> es(F_vbs_eigen);
  //
  //      std::cout << es.eigenvalues() << std::endl;
  //      auto ens_vbs_vir = es.eigenvalues().bottomRows(n_vbs - occ);
  //      ens = Eigen::VectorXd(n_vbs - n_frozen_core);
  //      ens << ens_occ, ens_vbs_vir;
  //
  //      std::cout << ens << std::endl;
  //      auto env = es.eigenvectors();
  //
  //      C_vir = env.rightCols(n_vbs - occ);
  //    }
  //
  //    // get all the trange1s
  //    auto tr_vir = tre->get_vir_tr1();
  //    auto tr_vbs = F_vbs.trange().data().back();
  //    detail::parallel_print_range_info(world, tr_vir, "Vir");
  //
  //    //    auto tr_all = tre->get_all_tr1();
  //    //    detail::parallel_print_range_info(world, tr_all, "Vbs");
  //
  //    auto C_vir_ta =
  //        array_ops::eigen_to_array<Tile>(world, C_vir, tr_vbs, tr_vir);
  //    auto vir_space =
  //        OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"),
  //        C_vir_ta);
  //    orbital_registry.add(vir_space);
  //
  //    auto mo_time1 = mpqc::fenced_now(world);
  //    auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  //    utility::print_par(world, "ClosedShell Dual Basis MO Build Time: ",
  //    mo_time,
  //                       " S\n");
  //
  //    return tre;
  //  }
  //
};

/**
 * RI-DBRMP2 class, only overide the compute function in DBRMP2 from RMP2
 */
template <typename Tile, typename Policy>
class RIDBRMP2 : public DBRMP2<Tile,Policy>{
public:
  /**
   * KeyVal constructor
   *
   * keywords, takes all keywords from DBRMP2
   */

  RIDBRMP2(const KeyVal& kv) : DBRMP2<Tile,Policy>(kv) {};
  ~RIDBRMP2() = default;
  double compute_scf_correction() override;
private:
  double compute() override;

};

extern template class DBRMP2<TA::TensorD, TA::SparsePolicy>;
extern template class RIDBRMP2<TA::TensorD, TA::SparsePolicy>;

}  // end of namespace mbpt
}  // end of namespace mpqc

#include "dbmp2_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_H_
