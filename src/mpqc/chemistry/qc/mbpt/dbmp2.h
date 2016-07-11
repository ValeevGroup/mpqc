//
// Created by Chong Peng on 6/13/16.
//

#ifndef MPQC_DBMP2_H
#define MPQC_DBMP2_H

#include <mpqc/chemistry/qc/mbpt/mp2.h>

namespace mpqc {
namespace mbpt {

template <typename Tile, typename Policy>
class DBMP2 {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegral = integrals::MolecularIntegral<Tile, Policy>;

  DBMP2(MolecularIntegral &mo_int, const rapidjson::Document &in) {
    auto mol = mo_int.atomic_integral().molecule();
    int occ = mol.occupation(0) / 2;
    Eigen::VectorXd orbital_energy;
    auto trange1_engine =
            closed_shell_dualbasis_mo_build_svd(mo_int, orbital_energy, in, mol, occ);

    mp2_ = std::make_shared<MP2<Tile, Policy>>(
        MP2<Tile, Policy>(mo_int, orbital_energy, trange1_engine));
  }

  double compute(const rapidjson::Document &in) {
    double mp2 = mp2_->compute(in);
    std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";
    int occ = mp2_->trange1_engine()->get_actual_occ();
    TArray F_ia;

    if(method == "four center"){
      F_ia = mp2_->mo_integral().compute(L"<i|F|a>");
    }
    else if(method == "df"){
      F_ia = mp2_->mo_integral().compute(L"<i|F|a>[df]");
    }
    else{
      throw std::runtime_error("Wrong MP2 Method");
    }

    double scf_correction = F_ia("i,a").reduce(ScfCorrection(mp2_->orbital_energy(),occ));

    if (F_ia.get_world().rank() == 0) {
      std::cout << "SCF Correction: " << scf_correction << std::endl;
    }
    return scf_correction + mp2;
  }

  private:

  struct ScfCorrection {
    using result_type = double;
    using argument_type = Tile;

    std::shared_ptr<Eig::VectorXd> vec_;
    unsigned int n_occ_;

    ScfCorrection(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
            : vec_(std::move(vec)), n_occ_(n_occ) { }

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
          me += (2.0 * data * data) / (e_ia);
        }
      }
    }
  };

  private:
//  std::shared_ptr<TRange1Engine> pulay_build(
//      integrals::MolecularIntegral<Tile, Policy> &mo_int,
//      OrbitalSpaceRegistry<TArray>& orbital_registry,
//      Eigen::VectorXd &ens, const rapidjson::Document &in,
//      const molecule::Molecule &mols, int occ) {
//    auto &ao_int = mo_int.atomic_integral();
//    auto &world = ao_int.world();
//    using TArray = TA::DistArray<Tile, Policy>;
//
//    auto mo_time0 = mpqc_time::fenced_now(world);
//    utility::print_par(world, "\nBuilding ClosedShell Dual Basis MO Orbital\n");
//
//    // solving occupied orbitals
//    TArray F;
//    if (ao_int.registry().have(Formula(L"<μ|F|ν>"))) {
//      F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>"));
//    } else {
//      F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
//    }
//
//    auto S = ao_int.compute(L"<κ|λ>");
//
//    MatrixD F_eig = array_ops::array_to_eigen(F);
//    MatrixD S_eig = array_ops::array_to_eigen(S);
//
//    // solve mo coefficients
//    Eig::GeneralizedSelfAdjointEigenSolver<MatrixD> es(F_eig, S_eig);
//
//    bool frozen_core =
//        in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
//    std::size_t n_frozen_core = 0;
//    if (frozen_core) {
//      n_frozen_core = mols.core_electrons();
//      utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons",
//                         "\n");
//      n_frozen_core = n_frozen_core / 2;
//    }
//
//    Eigen::VectorXd ens_occ =
//        es.eigenvalues().segment(n_frozen_core, occ - n_frozen_core);
//    std::cout << es.eigenvalues() << std::endl;
//    MatrixD C_all = es.eigenvectors();
//    MatrixD C_occ = C_all.block(0, 0, S_eig.rows(), occ);
//    MatrixD C_corr_occ =
//        C_all.block(0, n_frozen_core, S_eig.rows(), occ - n_frozen_core);
//
//    // finished solving occupied orbitals
//
//    // get all the sizes
//    auto n_obs = S.trange().elements().extent()[0];
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
//                                               vir_blocksize, n_frozen_core);
//    auto tr_obs = S.trange().data().back();
//    auto tr_occ = tre->compute_range(occ, occ_blocksize);
//    auto tr_corr_occ = tre->get_occ_tr1();
//    utility::parallel_print_range_info(world, tr_occ, "Occ");
//    utility::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
//
//    // convert to TA
//    auto C_occ_ta =
//        array_ops::eigen_to_array<Tile>(world, C_occ, tr_obs, tr_occ);
//    auto C_corr_occ_ta =
//        array_ops::eigen_to_array<Tile>(world, C_corr_occ, tr_obs, tr_corr_occ);
//
//    // insert to registry
//    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
//    auto occ_space =
//        OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
//    orbital_registry.add(occ_space);
//
//    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"),
//                                             OrbitalIndex(L"κ"), C_corr_occ_ta);
//    orbital_registry.add(corr_occ_space);
//
//    // solving the virtual orbitals
//    // find fock matrix
//    TArray F_vbs;
//    // if use density fitting
//    if (ao_int.registry().have(Formula(L"<μ|F|ν>[df]"))) {
//      F_vbs = ao_int.compute(Formula(L"<Α|F|Β>[df]"));
//    } else {
//      F_vbs = ao_int.compute(Formula(L"<Α|F|Β>"));
//    }
//    auto n_vbs = F_vbs.trange().elements().extent()[0];
//    tre = std::make_shared<TRange1Engine>(occ, n_vbs, occ_blocksize,
//                                          vir_blocksize, n_frozen_core);
//
//    // project Fock matrix against occ orbitals
//    MatrixD C_vir;
//    {
//      auto tr_vbs = F_vbs.trange().data().back();
//
//      // get identity matrix
//      Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(n_vbs, n_vbs);
//      auto identity_ta =
//          array_ops::eigen_to_array<Tile>(world, identity, tr_vbs, tr_vbs);
//      // get overlap
//      auto S_vbs = ao_int.compute(L"<Α|Β>");
//      // get density
//      MatrixD D_obs = C_occ * C_occ.transpose();
//
//      MatrixD D_vbs = MatrixD::Zero(n_vbs, n_vbs);
//      D_vbs.block(0, 0, n_obs, n_obs) << D_obs;
//      TArray D_vbs_ta =
//          array_ops::eigen_to_array<Tile>(world, D_vbs, tr_vbs, tr_vbs);
//
//      F_vbs("A,B") =
//          (identity_ta("A,D") - 0.5 * S_vbs("A,mu") * D_vbs_ta("mu,D")) *
//          F_vbs("D,E") *
//          (identity_ta("E,B") - 0.5 * D_vbs_ta("E,mu") * S_vbs("mu,B"));
//
//      MatrixD F_vbs_eigen = array_ops::array_to_eigen(F_vbs);
//
//      Eigen::SelfAdjointEigenSolver<MatrixD> es(F_vbs_eigen);
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
//    utility::parallel_print_range_info(world, tr_vir, "Vir");
//
//    //    auto tr_all = tre->get_all_tr1();
//    //    utility::parallel_print_range_info(world, tr_all, "Vbs");
//
//    auto C_vir_ta =
//        array_ops::eigen_to_array<Tile>(world, C_vir, tr_vbs, tr_vir);
//    auto vir_space =
//        OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
//    orbital_registry.add(vir_space);
//
//    auto mo_time1 = mpqc_time::fenced_now(world);
//    auto mo_time = mpqc_time::duration_in_s(mo_time0, mo_time1);
//    utility::print_par(world, "ClosedShell Dual Basis MO Build Time: ", mo_time,
//                       " S\n");
//
//    return tre;
//  }
//
//  std::shared_ptr<TRange1Engine> steele_build(
//          integrals::MolecularIntegral<Tile, Policy> &mo_int,
//      OrbitalSpaceRegistry<TArray>& orbital_registry,
//      Eigen::VectorXd &ens, const rapidjson::Document &in,
//      const molecule::Molecule &mols, int occ) {
//    auto &ao_int = mo_int.atomic_integral();
//    auto &world = ao_int.world();
//    using TArray = TA::DistArray<Tile, Policy>;
//
//    auto mo_time0 = mpqc_time::fenced_now(world);
//    utility::print_par(world, "\nBuilding ClosedShell Dual Basis MO Orbital\n");
//
//    // solving occupied orbitals
//    TArray F;
//    if (ao_int.registry().have(Formula(L"<μ|F|ν>"))) {
//      F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>"));
//    } else {
//      F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
//    }
//
//    auto S = ao_int.compute(L"<κ|λ>");
//
//    MatrixD F_eig = array_ops::array_to_eigen(F);
//    MatrixD S_eig = array_ops::array_to_eigen(S);
//
//    // solve mo coefficients
//    Eig::GeneralizedSelfAdjointEigenSolver<MatrixD> es(F_eig, S_eig);
//
//    bool frozen_core =
//        in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
//    std::size_t n_frozen_core = 0;
//    if (frozen_core) {
//      n_frozen_core = mols.core_electrons();
//      utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons",
//                         "\n");
//      n_frozen_core = n_frozen_core / 2;
//    }
//
//    std::cout << es.eigenvalues() << std::endl;
//    MatrixD C_occ = es.eigenvectors().leftCols(occ);
//
//    // finished solving occupied orbitals
//
//    // get all the sizes
//    auto n_obs = S.trange().elements().extent()[0];
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
//                                               vir_blocksize, n_frozen_core);
//    auto tr_obs = S.trange().data().back();
//    auto tr_occ = tre->compute_range(occ, occ_blocksize);
//    auto tr_corr_occ = tre->get_occ_tr1();
//    utility::parallel_print_range_info(world, tr_occ, "Occ");
//
//    // convert to TA
//    auto C_occ_ta =
//        array_ops::eigen_to_array<Tile>(world, C_occ, tr_obs, tr_occ);
//
//    // project to large basis set
//
//    TArray S_vbs_inv = ao_int.compute(L"<Α|Β>[inv]");
//    TArray S_vbs_obs = ao_int.compute(L"<Α|μ>");
//    auto n_vbs = S_vbs_inv.trange().elements().extent()[0];
//    auto tr_vbs = S_vbs_inv.trange().data().back();
//
//    TArray t;
//    t("A,mu") = S_vbs_inv("A,B") * S_vbs_obs("B,mu");
//
//
//    C_occ_ta("A,i") =  t("A,mu")* C_occ_ta("mu,i");
//
//    // insert to registry
//    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
//    auto occ_space =
//        OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"Α"), C_occ_ta);
//    orbital_registry.add(occ_space);
//
//    // solving the new orbitals in large basis set
//    // find fock matrix
//    TArray F_vbs;
//    // if use density fitting
//    if (ao_int.registry().have(Formula(L"<μ|F|ν>[df]"))) {
//      F_vbs = ao_int.compute(Formula(L"<Α|F|Β>[df]"));
//    } else {
//      F_vbs = ao_int.compute(Formula(L"<Α|F|Β>"));
//    }
//    tre = std::make_shared<TRange1Engine>(occ, n_vbs, occ_blocksize,
//                                          vir_blocksize, n_frozen_core);
//
//    MatrixD C_vir;
//    MatrixD C_corr_occ;
//    {
//      MatrixD F_vbs_eigen = array_ops::array_to_eigen(F_vbs);
//
//      Eigen::SelfAdjointEigenSolver<MatrixD> es(F_vbs_eigen);
//
//      assert(es.info() == Eigen::ComputationInfo::Success);
//      std::cout << es.eigenvalues() << std::endl;
//      ens = es.eigenvalues().bottomRows(n_vbs - n_frozen_core);
//      auto env = es.eigenvectors();
//
//      C_occ = env.leftCols(occ);
//      C_corr_occ = C_occ.rightCols(occ - n_frozen_core);
//      C_vir = env.rightCols(n_vbs - occ);
//    }
//
//    // get all the trange1s
//    auto tr_vir = tre->get_vir_tr1();
//    utility::parallel_print_range_info(world, tr_vir, "Vir");
//
//    //    auto tr_all = tre->get_all_tr1();
//    //    utility::parallel_print_range_info(world, tr_all, "Vbs");
//
//    C_occ_ta = array_ops::eigen_to_array<Tile>(world, C_occ, tr_vbs, tr_occ);
//
//    auto C_corr_occ_ta =
//        array_ops::eigen_to_array<Tile>(world, C_corr_occ, tr_vbs, tr_corr_occ);
//
//    auto C_vir_ta =
//        array_ops::eigen_to_array<Tile>(world, C_vir, tr_vbs, tr_vir);
//
//    // insert to registry
//    occ_space =
//        OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"Α"), C_occ_ta);
//    orbital_registry.remove(OrbitalIndex(L"m"));
//    orbital_registry.add(occ_space);
//
//    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"),
//                                             OrbitalIndex(L"Α"), C_corr_occ_ta);
//    orbital_registry.add(corr_occ_space);
//
//    auto vir_space =
//        OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
//    orbital_registry.add(vir_space);
//
//    auto mo_time1 = mpqc_time::fenced_now(world);
//    auto mo_time = mpqc_time::duration_in_s(mo_time0, mo_time1);
//    utility::print_par(world, "ClosedShell Dual Basis MO Build Time: ", mo_time,
//                       " S\n");
//
//    return tre;
//  }


  std::shared_ptr<TRange1Engine> closed_shell_dualbasis_mo_build_svd(
          integrals::MolecularIntegral<Tile, Policy> &mo_int,
          Eigen::VectorXd &ens, const rapidjson::Document &in,
          const molecule::Molecule &mols, int occ) {
    auto &ao_int = mo_int.atomic_integral();
    auto &world = ao_int.world();
    using TArray = TA::DistArray<Tile, Policy>;

    auto mo_time0 = mpqc_time::fenced_now(world);
    utility::print_par(world, "\nBuilding ClosedShell Dual Basis MO Orbital\n");

    // solving occupied orbitals
    TArray F;
    if (ao_int.registry().have(Formula(L"<μ|F|ν>"))) {
      F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>"));
    } else {
      F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
    }

    auto S = ao_int.compute(L"<κ|λ>");

    MatrixD F_eig = array_ops::array_to_eigen(F);
    MatrixD S_eig = array_ops::array_to_eigen(S);

    // solve mo coefficients
    Eig::GeneralizedSelfAdjointEigenSolver<MatrixD> es(F_eig, S_eig);

    bool frozen_core =
            in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
    std::size_t n_frozen_core = 0;
    if (frozen_core) {
      n_frozen_core = mols.core_electrons();
      utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons",
                         "\n");
      n_frozen_core = n_frozen_core / 2;
    }

    Eigen::VectorXd ens_occ = es.eigenvalues().segment(n_frozen_core, occ - n_frozen_core);
    std::cout << ens_occ << std::endl;
    MatrixD C_all = es.eigenvectors();
    MatrixD C_occ = C_all.block(0, 0, S_eig.rows(), occ);
    MatrixD C_corr_occ = C_all.block(0, n_frozen_core, S_eig.rows(), occ - n_frozen_core);

    // finished solving occupied orbitals


    // start to solve virtual orbitals
    auto S_vbs = ao_int.compute(L"<Α|Β>");
    auto S_obs_vbs = ao_int.compute(L"<μ|Α>");

    // construct C_vbs
    MatrixD C_vbs;
    std::size_t nbf_vbs;
    std::size_t nbf_v;
    {

      // S_A^B -(1/2)
      MatrixD S_vbs_eigen = array_ops::array_to_eigen(S_vbs);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S_vbs_eigen);
      MatrixD X_vbs_eigen_inv = es2.operatorInverseSqrt();

      // orthogonalize
      MatrixD S_obs_vbs_eigen = array_ops::array_to_eigen(S_obs_vbs);
      MatrixD S_obs_vbs_ortho_eigen = S_obs_vbs_eigen * X_vbs_eigen_inv;

      MatrixD X_i_mu = C_occ.transpose()*S_obs_vbs_ortho_eigen;

      // SVD solve
      Eigen::JacobiSVD<MatrixD> svd(X_i_mu, Eigen::ComputeFullV);
      MatrixD V_eigen = svd.matrixV();
      nbf_vbs = S_obs_vbs_ortho_eigen.cols();
      nbf_v = nbf_vbs - svd.nonzeroSingularValues();
      std::cout << "Virtual Orbital Size:  " << nbf_v << std::endl;
      MatrixD Vnull(nbf_vbs, nbf_v);
      Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_vbs, nbf_v);
      C_vbs = X_vbs_eigen_inv.transpose() * Vnull;

    }


    // get all the sizes
    auto n_obs = S.trange().elements().extent()[0];
    std::size_t mo_blocksize = in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
    std::size_t vir_blocksize = in.HasMember("VirBlockSize")
                                ? in["VirBlockSize"].GetInt()
                                : mo_blocksize;
    std::size_t occ_blocksize = in.HasMember("OccBlockSize")
                                ? in["OccBlockSize"].GetInt()
                                : mo_blocksize;
    utility::print_par(world, "OccBlockSize: ", occ_blocksize, "\n");
    utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");

    auto tre = std::make_shared<TRange1Engine>(occ, nbf_vbs, occ_blocksize,
                                               vir_blocksize, n_frozen_core);
    auto tr_obs = S.trange().data().back();
    auto tr_vbs = S_vbs.trange().data().back();
    auto tr_occ = tre->compute_range(occ, occ_blocksize);
    auto tr_corr_occ = tre->get_occ_tr1();
    auto tr_vir = tre->get_vir_tr1();

    utility::parallel_print_range_info(world, tr_occ, "Occ");
    utility::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
    utility::parallel_print_range_info(world, tr_vir, "Vir");

    // convert to TA
    auto C_occ_ta = array_ops::eigen_to_array<Tile>(world, C_occ, tr_obs, tr_occ);
    auto C_corr_occ_ta = array_ops::eigen_to_array<Tile>(world, C_corr_occ, tr_obs, tr_corr_occ);
    auto C_vir_ta = array_ops::eigen_to_array<Tile>(world, C_vbs, tr_vbs, tr_vir);

    // insert to registry
    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
    auto occ_space = OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
    mo_int.orbital_space()->add(occ_space);

    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
    mo_int.orbital_space()->add(corr_occ_space);

    auto vir_space = OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
    mo_int.orbital_space()->add(vir_space);


    // solve energy in virtual orbital
    TArray F_vbs;
    if (ao_int.registry().have(Formula(L"<μ|F|ν>"))) {
      F_vbs = mo_int.compute(L"<a|F|b>");
    } else {
      F_vbs = mo_int.compute(L"<a|F|b>[df]");
    }
    MatrixD F_vbs_mo_eigen = array_ops::array_to_eigen(F_vbs);
//    std::cout << "F_vbs MO" << std::endl;
//    std::cout << F_vbs_mo_eigen << std::endl;

    Eigen::SelfAdjointEigenSolver<MatrixD> es3(F_vbs_mo_eigen);
    auto ens_vir = es3.eigenvalues();

    ens = Eigen::VectorXd(nbf_vbs - n_frozen_core);
    ens << ens_occ, ens_vir;

    std::cout << "Energy of Orbitals " << std::endl;
    std::cout << ens << std::endl;

    // resolve the virtual orbitals
    MatrixD C_vir_rotate = es3.eigenvectors();
    C_vbs = C_vbs*C_vir_rotate;
    TArray C_vir_ta_new = array_ops::eigen_to_array<Tile>(world, C_vbs, tr_vbs, tr_vir);

    // remove old virtual orbitals
    mo_int.orbital_space()->remove(OrbitalIndex(L"a"));
    mo_int.registry().remove_orbital(world,OrbitalIndex(L"a"));

    // add new virtual orbial
    vir_space = OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta_new);
    mo_int.orbital_space()->add(vir_space);

    return tre;
  }

 private:
  std::shared_ptr<MP2<Tile, Policy>> mp2_;
};

}  // end of namespace mbpt
}  // end of namespace mpqc

#endif  // MPQC_DBMP2_H
