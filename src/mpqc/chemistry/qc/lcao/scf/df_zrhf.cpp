#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"

#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"

namespace mpqc {
namespace lcao {

DFzRHF::DFzRHF(const KeyVal &kv) : zRHF(kv) {}

void DFzRHF::init_other() {
  auto &ao_factory = this->ao_factory();
  auto &world = ao_factory.world();

  auto t0_init = mpqc::fenced_now(world);
  mpqc::time_point t0, t1;

  // overlap matrix
  M_ = this->S_;

  // 1*N charge matrix of auxiliary basis within one unit cell
  auto unit_basis = gaussian::detail::create_unit_basis();
  ao_factory.basis_registry()->add(OrbitalIndex(L"U"), unit_basis);
  auto charge_mat = ao_factory.compute(L"< U | Κ >");
  // charge vector of auxiliary basis within one unit cell
  auto charge_vec = gaussian::detail::take_row_from_2D_array(charge_mat, 0);
  double charge_2 = charge_vec("X") * charge_vec("X");
  q_ = std::sqrt(charge_2);
  // normalized charge vector of auxiliary basis within one unit cell
  n_("X") = (1.0 / q_) * charge_vec("X");

  // projection matrix that projects matrix X on to charge vector
  P_para_("X, Y") = n_("X") * n_("Y");

  // projection matrix that projects matrix X orthogonal to charge vector
  identity_ = ao_factory.compute(L"< Κ | I | Λ >");
  P_perp_("X, Y") = identity_("X, Y") - P_para_("X, Y");

  // 2-body 2-center integrals
  V_ = ao_factory.compute(L"( Κ | G | Λ )");

  // perpendicular part of 2-body 2-center integrals
  t0 = mpqc::fenced_now(world);
  V_perp_("X, Y") = P_perp_("X, Z") * V_("Z, W") * P_perp_("W, Y");
  t1 = mpqc::fenced_now(world);
  auto t_v_perp = mpqc::duration_in_s(t0, t1);

  // A inverse where A = V_perp + P_para
  t0 = mpqc::fenced_now(world);
  TArray A;
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
  TArray L_inv =
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

  Gamma_vec_ = ao_factory.compute_direct_vector(L"( Κ | G|κ λ)");

  auto t1_init = mpqc::fenced_now(world);
  double t_init_other = mpqc::duration_in_s(t0_init, t1_init);

  if (this->print_detail_) {
    ExEnv::out0() << "\nRI-J init time decomposition:\n"
                  << "\tV perp:              " << t_v_perp << " s\n"
                  << "\tA = V_perp + P_para: " << t_a << " s\n"
                  << "\tL inv:               " << t_l_inv << " s\n"
                  << "\tA inv:               " << t_a_inv << " s\n"
                  << "\tIM:                  " << t_im << " s" << std::endl;
  }
  ExEnv::out0() << "\nInit RI-J time:      " << t_init_other << " s\n"
                << std::endl;
}

DFzRHF::TArray DFzRHF::J_builder() {
  auto &ao_factory = this->ao_factory();
  auto &world = ao_factory.world();

  mpqc::time_point t0, t1;
  auto t0_j_builder = mpqc::fenced_now(world);

  // 3-center 2-electron direct integrals contracted with density matrix
  G_ = ao_factory.compute_direct(L"( Κ | G|κ λ)");

  // Build [CD]_X = C_Xμν D_μν
  double t_w_para, t_w;
  {
    // intermediate for C_para_Xμν D_μν
    t0 = mpqc::fenced_now(world);
    TArray interm;
    interm("X") = (1.0 / q_) * M_("mu, nu") * n_("X") * this->D_("mu, nu");
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
  TArray J, J_part1, J_part2;
  TArray VCD;  // just an intermediate
  VCD("X") = V_("X, Y") * CD_("Y");

  // Build Part #1 of J_μν
  double t_j1_interm, t_j1_contr, t_j1;
  {
    t0 = mpqc::fenced_now(world);

    auto t0_j1_interm = mpqc::fenced_now(world);
    TArray interm;
    interm("X") = 2.0 * CD_("X") - IP_("Y, X") * VCD("Y");
    auto t1_j1_interm = mpqc::fenced_now(world);
    t_j1_interm = mpqc::duration_in_s(t0_j1_interm, t1_j1_interm);

    auto t0_j1_contr = mpqc::fenced_now(world);
    auto RJ_size = ao_factory.RJ_size();
    for (auto RJ = 0; RJ < RJ_size; ++RJ) {
      auto &g = Gamma_vec_[RJ];
      if (RJ == 0)
        J_part1("mu, nu") = g("X, mu, nu") * interm("X");
      else
        J_part1("mu, nu") += g("X, mu, nu") * interm("X");
    }
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
    TArray interm;
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

}  // namespace lcao

}  // namespace mpqc
