#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"

#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"

namespace mpqc {
namespace lcao {

DFzRHF::DFzRHF(const KeyVal &kv) : zRHF(kv) {}

void DFzRHF::init_other() {
    auto &ao_factory = this->ao_factory();
    auto &world = ao_factory.world();

    if (this->print_detail_) {
        ExEnv::out0() << "\nRI-J initialization:" << std::endl;
    }

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
    TArray identity = ao_factory.compute(L"< Κ | I | Λ >");
    P_perp_("X, Y") = identity("X, Y") - P_para_("X, Y");

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

    Gamma_vec_ = ao_factory.compute_direct_vector(L"( Κ | G|κ λ)");

    auto t1_init = mpqc::fenced_now(world);
    double t_init_other = mpqc::duration_in_s(t0_init, t1_init);

    if (this->print_detail_) {
        ExEnv::out0() << "\tV perp:              " << t_v_perp << " s" << std::endl;
        ExEnv::out0() << "\tA = V_perp + P_para: " << t_a << " s" << std::endl;
        ExEnv::out0() << "\tL inv:               " << t_l_inv << " s" << std::endl;
        ExEnv::out0() << "\tA inv:               " << t_a_inv << " s" << std::endl;
        ExEnv::out0() << "\nInit RI-J time:      " << t_init_other << " s\n" << std::endl;
    }
}

DFzRHF::TArray DFzRHF::J_builder() {
    auto &ao_factory = this->ao_factory();
    auto &world = ao_factory.world();

    auto t0_j_builder = mpqc::fenced_now(world);

    // G_X: 2-body 3-center direct integrals contracted with density matrix
    G_ = ao_factory.compute(L"( Κ | G|κ λ)");

    // intermediate for C_para_Xμν D_μν
    W_para_("X") = (1.0 / q_) * M_("mu, nu") * n_("X") * this->D_("mu, nu");

    // intermediate for C_Xμν D_μν
    W_("X") = inv_("X, Y") * P_perp_("Y, Z") * (G_("Z") - V_("Z, W") * W_para_("W"));

    // build DF Coulomb term
    TArray J, J_part1, J_part2;
    t0 = mpqc::fenced_now(world);
    auto RJ_size = ao_factory.RJ_size();
    for (auto RJ = 0; RJ < RJ_size; ++RJ) {
        auto &g = Gamma_vec_[RJ];
        if (RJ == 0)
            J_part1("mu, nu") = g("X, mu, nu") * W_("X");
        else
            J_part1("mu, nu") += g("X, mu, nu") * W_("X");
    }
    t1 = mpqc::fenced_now(world);
    auto t_j1 = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    TArray VW;
    VW("X") = V_("X, Y") * W_("Y");
    J_part2("mu, nu") = C_df_("X, mu, nu") * VW("X");
    t1 = mpqc::fenced_now(world);
    auto t_j2 = mpqc::duration_in_s(t0, t1);




    mpqc::time_point t0, t1;
    // DF coeffs parallel to charge vector
    t0 = mpqc::fenced_now(world);
    C_para_("X, mu, nu") =
        (1.0 / q_) * M_("mu, nu") * n_("X");
    t1 = mpqc::fenced_now(world);
    auto t_c_para = mpqc::duration_in_s(t0, t1);

    // perpendicular part of 3-body 2-center integrals
    t0 = mpqc::fenced_now(world);
    TArray Gamma_perp;
    Gamma_perp("X, mu, nu") = P_perp_("X, Y") * Gamma_("Y, mu, nu");
    t1 = mpqc::fenced_now(world);
    auto t_gamma_perp = mpqc::duration_in_s(t0, t1);


    // perpendicular η
    t0 = mpqc::fenced_now(world);
    TArray Eta_perp;
    Eta_perp("X, mu, nu") = P_perp_("X, Y") * V_("Y, Z") * C_para_("Z, mu, nu");
    t1 = mpqc::fenced_now(world);
    auto t_eta_perp = mpqc::duration_in_s(t0, t1);


    // DF coeffs orthogonal to charge vector
    t0 = mpqc::fenced_now(world);
    C_perp_("X, mu, nu") =
        inv_("X, Y") * (Gamma_perp("Y, mu, nu") - Eta_perp("Y, mu, nu"));
    t1 = mpqc::fenced_now(world);
    auto t_c_perp = mpqc::duration_in_s(t0, t1);

    // total DF coeffs
    t0 = mpqc::fenced_now(world);
    C_df_("X, mu, nu") = C_para_("X, mu, nu") + C_perp_("X, mu, nu");
    t1 = mpqc::fenced_now(world);
    auto t_c = mpqc::duration_in_s(t0, t1);

    // build DF Coulomb term
    TArray J, J_part1, J_part2, CD;
    t0 = mpqc::fenced_now(world);
    CD("X") = C_df_("X, rho, sig") * this->D_("rho, sig");
    t1 = mpqc::fenced_now(world);
    auto t_cd = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    J_part1("mu, nu") = 2.0 * Gamma_("X, mu, nu") * CD("X");
    t1 = mpqc::fenced_now(world);
    auto t_j1 = mpqc::duration_in_s(t0, t1);

    t0 = mpqc::fenced_now(world);
    TArray VCD;
    VCD("X") = V_("X, Y") * CD("Y");
    J_part2("mu, nu") = C_df_("X, mu, nu") * VCD("X");
    t1 = mpqc::fenced_now(world);
    auto t_j2 = mpqc::duration_in_s(t0, t1);

//    J_part1("mu, nu") = 2.0 * Gamma_("X, mu, nu") * C_df_("X, rho, sig") * this->D_("rho, sig");
//    J_part2("mu, nu") = C_df_("X, mu, nu") * V_("X, Y") * C_df_("Y, rho, sig") * this->D_("rho, sig");
    J("mu, nu") = J_part1("mu, nu") - J_part2("mu, nu");

    auto t1_j_builder = mpqc::fenced_now(world);
    auto t_tot = mpqc::duration_in_s(t0_j_builder, t1_j_builder);

    if (this->print_detail_) {
        ExEnv::out0() << "\nRI-J timing decomposition:" << std::endl;
        ExEnv::out0() << "\tC para:              " << t_c_para << " s" << std::endl;
        ExEnv::out0() << "\tGamma perp:          " << t_gamma_perp << " s" << std::endl;
        ExEnv::out0() << "\tEta perp:            " << t_eta_perp << " s" << std::endl;
        ExEnv::out0() << "\tC perp:              " << t_c_perp << " s" << std::endl;
        ExEnv::out0() << "\tC = C_para + C_perp: " << t_c << " s" << std::endl;
        ExEnv::out0() << "\tCD = C_Xμν D_μν:     " << t_cd << " s" << std::endl;
        ExEnv::out0() << "\tJ_part1:             " << t_j1 << " s" << std::endl;
        ExEnv::out0() << "\tJ_part2:             " << t_j2 << " s" << std::endl;
        ExEnv::out0() << "\nTotal J builder time: " << t_tot << " s" << std::endl;
    }
    return J;
}

}  // namespace lcao

}  // namespace mpqc
