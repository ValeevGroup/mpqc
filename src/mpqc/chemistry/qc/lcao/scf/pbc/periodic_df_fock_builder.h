#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_DF_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_DF_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_three_center_contraction_builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicDFFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using DirectTArray = typename Factory::DirectTArray;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;
  using PFC_Builder = PeriodicFourCenterFockBuilder<Tile, Policy>;

  PeriodicDFFockBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    print_detail_ = ao_factory_.print_detail();
    auto &world = ao_factory_.world();

    auto t0_init = mpqc::fenced_now(world);
    mpqc::time_point t0, t1;

    // overlap matrix
    M_ = ao_factory_.compute(L"<κ|λ>");

    // 1*N charge matrix of auxiliary basis within one unit cell
    auto unit_basis = ::mpqc::lcao::gaussian::detail::create_unit_basis();
    ao_factory_.basis_registry()->add(OrbitalIndex(L"U"), unit_basis);
    auto charge_mat = ao_factory_.compute(L"< U | Κ >");
    // charge vector of auxiliary basis within one unit cell
    auto charge_vec =
        ::mpqc::lcao::gaussian::detail::take_row_from_2D_array(charge_mat, 0);
    double charge_2 = charge_vec("X") * charge_vec("X");
    q_ = std::sqrt(charge_2);
    // normalized charge vector of auxiliary basis within one unit cell
    n_("X") = (1.0 / q_) * charge_vec("X");

    // projection matrix that projects matrix X on to charge vector
    P_para_("X, Y") = n_("X") * n_("Y");

    // projection matrix that projects matrix X orthogonal to charge vector
    identity_ = ao_factory_.compute(L"< Κ | I | Λ >");
    P_perp_("X, Y") = identity_("X, Y") - P_para_("X, Y");

    // 2-body 2-center integrals
    V_ = ao_factory_.compute(L"( Κ | G | Λ )");

    // perpendicular part of 2-body 2-center integrals
    t0 = mpqc::fenced_now(world);
    V_perp_("X, Y") = P_perp_("X, Z") * V_("Z, W") * P_perp_("W, Y");
    t1 = mpqc::fenced_now(world);
    auto t_v_perp = mpqc::duration_in_s(t0, t1);

    // A inverse where A = V_perp + P_para
    t0 = mpqc::fenced_now(world);
    array_type A;
    A("X, Y") = V_perp_("X, Y") + P_para_("X, Y");
    t1 = mpqc::fenced_now(world);
    auto t_a = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    auto A_eig = array_ops::array_to_eigen(A);
    using MatType = decltype(A_eig);
    MatType L_inv_eig = MatType(Eigen::LLT<MatType>(A_eig).matrixL()).inverse();
    auto tr0 = A.trange().data()[0];
    auto tr1 = A.trange().data()[1];
    assert(tr0 == tr1 && "Matrix A = LLT must be symmetric!");
    array_type L_inv =
        array_ops::eigen_to_array<Tile, Policy>(world, L_inv_eig, tr0, tr1);
    t1 = mpqc::fenced_now(world);
    auto t_l_inv = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    inv_("X, Y") = L_inv("Z, X") * L_inv("Z, Y");
    t1 = mpqc::fenced_now(world);
    auto t_a_inv = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    IP_("X, Y") = inv_("X, Z") * P_perp_("Z, Y");
    t1 = mpqc::fenced_now(world);
    auto t_im = mpqc::duration_in_s(t0, t1);

    // collect information to construct 3-center and 4-center builders
    auto basis = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    auto aux_basis = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));
    auto dcell = ao_factory_.unitcell().dcell();
    auto R_max = ao_factory_.R_max();
    auto RJ_max = ao_factory_.RJ_max();
    auto RD_max = ao_factory_.RD_max();
    auto R_size = ao_factory_.R_size();
    auto RJ_size = ao_factory_.RJ_size();
    auto RD_size = ao_factory_.RD_size();
    auto screen = ao_factory_.screen();
    auto screen_threshold = ao_factory_.screen_threshold();

    // construct PerioidcThreeCenterContractionBuilder for contractions
    // involving 3-center ints
    three_center_builder_ = std::make_unique<PTC_Builder>(
        world, basis, aux_basis, dcell, R_max, RJ_max, RD_max, R_size, RJ_size,
        RD_size, screen, screen_threshold);

    auto t1_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_init, t1_init);

    auto t0_k_init = mpqc::fenced_now(world);
    // construct PerioidcFourCenterFockBuilder for exchange term
    k_builder_ = std::make_unique<PFC_Builder>(
        world, basis, basis, dcell, R_max, RJ_max, RD_max, R_size, RJ_size,
        RD_size, false, true, screen, screen_threshold);
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    if (this->print_detail_) {
      ExEnv::out0() << "\nRI-J init time decomposition:\n"
                    << "\tV perp:              " << t_v_perp << " s\n"
                    << "\tA = V_perp + P_para: " << t_a << " s\n"
                    << "\tL inv:               " << t_l_inv << " s\n"
                    << "\tA inv:               " << t_a_inv << " s\n"
                    << "\tIM:                  " << t_im << " s" << std::endl;
    }
    ExEnv::out0() << "\nInit RI-J time:      " << t_j_init << " s" << std::endl;
    ExEnv::out0() << "\nInit Four-Center-K time:      " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicDFFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision) override {
    array_type G;

    // the '-' sign is embeded in K builder
    G("mu, nu") = 2.0 * compute_J(D, target_precision)("mu, nu") +
                  compute_K(D, target_precision)("mu, nu");

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

 private:
  Factory &ao_factory_;
  bool print_detail_;
  std::unique_ptr<PTC_Builder> three_center_builder_;
  std::unique_ptr<PFC_Builder> k_builder_;

  array_type M_;         // charge matrix of product density <μ|ν>
  array_type n_;         // normalized charge vector <Κ>
  double q_;             // total charge of auxiliary basis functions
  array_type P_para_;    // projection matrix that projects X onto auxiliary
                         // charge vector
  array_type P_perp_;    // projection matrix that projects X onto the subspace
                         // orthogonal to auxiliary charge vector
  array_type V_;         // 2-center 2-electron integrals
  array_type V_perp_;    // part of 2-center 2-electron integrals that is
                         // orthogonal to auxiliary charge vector
  array_type G_;         // 3-center 2-electron direct integrals contracted with
                         // density matrix
  array_type inv_;       // A inverse where A = V_perp + P_para
  array_type identity_;  // idensity matrix
  std::vector<DirectTArray> Gamma_vec_;  // vector of 3-center 2-electron direct
                                         // integrals. vector size = RJ_size_
  array_type CD_;                        // intermediate for C_Xμν D_μν
  array_type IP_;                        // intermediate for inv_XY P_perp_YZ

 private:
  array_type compute_J(const array_type &D, double target_precision) {
    auto &world = ao_factory_.world();
    // feed density matrix to Factory
    ao_factory_.set_density(D);

    mpqc::time_point t0, t1;
    auto t0_j_builder = mpqc::fenced_now(world);

    // 3-center 2-electron direct integrals contracted with density matrix
    t0 = mpqc::fenced_now(world);
    G_ = three_center_builder_->template contract_with<1>(D, target_precision);
    t1 = mpqc::fenced_now(world);
    auto t_3c_d_contr = mpqc::duration_in_s(t0, t1);

    // Build [CD]_X = C_Xμν D_μν
    double t_w_para, t_w;
    {
      // intermediate for C_para_Xμν D_μν
      t0 = mpqc::fenced_now(world);
      array_type interm;
      interm("X") = (1.0 / q_) * M_("mu, nu") * n_("X") * D("mu, nu");
      t1 = mpqc::fenced_now(world);
      t_w_para = mpqc::duration_in_s(t0, t1);

      // intermediate for C_Xμν D_μν
      t0 = mpqc::fenced_now(world);
      CD_ = interm;
      CD_("X") += IP_("X, Y") * (G_("Y") - V_("Y, Z") * interm("Z"));
      t1 = mpqc::fenced_now(world);
      t_w = mpqc::duration_in_s(t0, t1);
    }

    // Build DF Coulomb term J_μν
    array_type J, J_part1, J_part2;
    array_type VCD;  // just an intermediate
    VCD("X") = V_("X, Y") * CD_("Y");

    // Build Part #1 of J_μν
    double t_j1_interm, t_j1_contr, t_j1;
    {
      t0 = mpqc::fenced_now(world);

      auto t0_j1_interm = mpqc::fenced_now(world);
      array_type interm;
      interm("X") = 2.0 * CD_("X") - IP_("Y, X") * VCD("Y");
      auto t1_j1_interm = mpqc::fenced_now(world);
      t_j1_interm = mpqc::duration_in_s(t0_j1_interm, t1_j1_interm);

      auto t0_j1_contr = mpqc::fenced_now(world);
      J_part1 = three_center_builder_->template contract_with<2>(
          interm, target_precision);
      auto t1_j1_contr = mpqc::fenced_now(world);
      t_j1_contr = mpqc::duration_in_s(t0_j1_contr, t1_j1_contr);

      t1 = mpqc::fenced_now(world);
      t_j1 = mpqc::duration_in_s(t0, t1);
    }

    // Build Part #2 of J_μν
    double t_j2_interm, t_j2_eq, t_j2;
    {
      t0 = mpqc::fenced_now(world);

      auto t0_j2_interm = mpqc::fenced_now(world);
      array_type interm;
      interm("X, Y") = IP_("X, Z") * V_("Z, Y") - identity_("X, Y");
      auto t1_j2_interm = mpqc::fenced_now(world);
      t_j2_interm = mpqc::duration_in_s(t0_j2_interm, t1_j2_interm);

      auto t0_j2_eq = mpqc::fenced_now(world);
      double prefactor = interm("X, Z") * n_("Z") * VCD("X");
      J_part2("mu, nu") = (prefactor / q_) * M_("mu, nu");
      auto t1_j2_eq = mpqc::fenced_now(world);
      t_j2_eq = mpqc::duration_in_s(t0_j2_eq, t1_j2_eq);

      t1 = mpqc::fenced_now(world);
      t_j2 = mpqc::duration_in_s(t0, t1);
    }

    // Build J_μν
    J("mu, nu") = J_part1("mu, nu") + J_part2("mu, nu");

    auto t1_j_builder = mpqc::fenced_now(world);
    auto t_tot = mpqc::duration_in_s(t0_j_builder, t1_j_builder);

    if (this->print_detail_) {
      ExEnv::out0() << "\nRI-J timing decomposition:\n"
                    << "\tSum_RJ (X|μν) D_μν:   " << t_3c_d_contr << " s\n"
                    << "\tC_para_Xμν D_μν:      " << t_w_para << " s\n"
                    << "\tC_Xμν D_μν:           " << t_w << " s\n"
                    << "\tJ_part1:              " << t_j1 << " s\n"
                    << "\t  interm:             " << t_j1_interm << " s\n"
                    << "\t  contr:              " << t_j1_contr << " s\n"
                    << "\tJ_part2:              " << t_j2 << " s\n"
                    << "\t  interm :            " << t_j2_interm << " s\n"
                    << "\t  eq:                 " << t_j2_eq << " s\n"
                    << "\nTotal J builder time: " << t_tot << " s" << std::endl;
    }
    return J;
  }

  array_type compute_K(const array_type &D, double target_precision) {
    return k_builder_->operator()(D, target_precision);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_DF_FOCK_BUILDER_H_
