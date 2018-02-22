#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_three_center_contraction_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"
#include "mpqc/math/linalg/cholesky_inverse.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicRIJBuilder {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;

  PeriodicRIJBuilder(Factory &ao_factory) : ao_factory_(ao_factory) { init(); }

  PeriodicRIJBuilder(Factory &ao_factory, const Vector3i &rj_max)
      : ao_factory_(ao_factory) {
    ao_factory_.set_rjmax(rj_max);
    init();
  }

  ~PeriodicRIJBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    return compute_J(D, target_precision);
  }

 private:
  Factory &ao_factory_;
  bool print_detail_;
  std::unique_ptr<PTC_Builder> three_center_builder_;

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
  array_type CD_;        // intermediate for C_Xμν D_μν
  array_type IP_;        // intermediate for inv_XY P_perp_YZ

 private:
  void init() {
    print_detail_ = ao_factory_.print_detail();
    auto &world = ao_factory_.world();

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
    inv_ = array_ops::eigen_inverse(A);
    t1 = mpqc::fenced_now(world);
    auto t_a_inv = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    IP_("X, Y") = inv_("X, Z") * P_perp_("Z, Y");
    t1 = mpqc::fenced_now(world);
    auto t_im = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    three_center_builder_ = std::make_unique<PTC_Builder>(ao_factory_);
    t1 = mpqc::fenced_now(world);
    auto t_3c_ctor = mpqc::duration_in_s(t0, t1);

    if (this->print_detail_) {
      ExEnv::out0() << "\nRI-J init time decomposition:\n"
                    << "\tV perp:              " << t_v_perp << " s\n"
                    << "\tA = V_perp + P_para: " << t_a << " s\n"
                    << "\tA inv:               " << t_a_inv << " s\n"
                    << "\tIM:                  " << t_im << " s\n"
                    << "\t3-c builder ctor:    " << t_3c_ctor << " s"
                    << std::endl;
    }
  }

  array_type compute_J(const array_type &D, double target_precision) {
    auto &world = ao_factory_.world();

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
      double prefactor;
      {
        const auto R_max = ao_factory_.R_max();
        const auto RD_max = ao_factory_.RD_max();
        prefactor = ::mpqc::pbc::detail::dot_product(M_, D, R_max, RD_max);
      }
      interm("X") = (prefactor / q_) * n_("X");
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
    double t_j1;
    {
      t0 = mpqc::fenced_now(world);
      array_type interm;
      interm("X") = 2.0 * CD_("X") - IP_("Y, X") * VCD("Y");

      J_part1 = three_center_builder_->template contract_with<2>(
          interm, target_precision);
      t1 = mpqc::fenced_now(world);
      t_j1 = mpqc::duration_in_s(t0, t1);
    }

    // Build Part #2 of J_μν
    double t_j2;
    {
      t0 = mpqc::fenced_now(world);
      array_type interm;
      interm("X, Y") = IP_("X, Z") * V_("Z, Y") - identity_("X, Y");

      double prefactor = interm("X, Z") * n_("Z") * VCD("X");
      J_part2("mu, nu") = (prefactor / q_) * M_("mu, nu");
      t1 = mpqc::fenced_now(world);
      t_j2 = mpqc::duration_in_s(t0, t1);
    }

    // Build J_μν
    J("mu, nu") = J_part1("mu, nu") + J_part2("mu, nu");

    auto t1_j_builder = mpqc::fenced_now(world);
    auto t_tot = mpqc::duration_in_s(t0_j_builder, t1_j_builder);

    if (this->print_detail_) {
      ExEnv::out0() << "\nRI-J time decomposition:\n"
                    << "\tSum_RJ (X|μν) D_μν:   " << t_3c_d_contr << " s\n"
                    << "\tC_para_Xμν D_μν:      " << t_w_para << " s\n"
                    << "\tC_Xμν D_μν:           " << t_w << " s\n"
                    << "\tJ_part1:              " << t_j1 << " s\n"
                    << "\tJ_part2:              " << t_j2 << " s\n"
                    << "\nTotal J builder time: " << t_tot << " s" << std::endl;
    }
    return J;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_BUILDER_H_
