#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"

#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"

namespace mpqc {
namespace lcao {

DFzRHF::DFzRHF(const KeyVal &kv) : zRHF(kv) {}

void DFzRHF::init_other() {
    auto &ao_factory = this->ao_factory();

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

    // 2-body 2-center integrals
    V_ = ao_factory.compute(L"( Κ | G | Λ )");

    // 2-body 3-center integrals
    // TODO use direct integral
    Gamma_ = ao_factory.compute(L"( Κ | G|κ λ)");

    // projection matrix that projects matrix X on to charge vector
    P_para_("X, Y") = n_("X") * n_("Y");

    // projection matrix that projects matrix X orthogonal to charge vector
    TArray identity = ao_factory.compute(L"< Κ | I | Λ >");
    P_perp_("X, Y") = identity("X, Y") - P_para_("X, Y");

}

DFzRHF::TArray DFzRHF::J_builder() {
    auto &ao_factory = this->ao_factory();
    auto &world = ao_factory.world();

    // DF coeffs parallel to charge vector
    C_para_("X, mu, nu") =
        (1.0 / q_) * M_("mu, nu") * n_("X");

    // perpendicular part of 3-body 2-center integrals
    TArray Gamma_perp;
    Gamma_perp("X, mu, nu") = P_perp_("X, Y") * Gamma_("Y, mu, nu");

    // perpendicular part of 2-body 2-center integrals
    TArray V_perp;
    V_perp("X, Y") = P_perp_("X, Z") * V_("Z, W") * P_perp_("W, Y");

    // perpendicular η
    TArray Eta_perp;
    Eta_perp("X, mu, nu") = P_perp_("X, Y") * V_("Y, Z") * C_para_("Z, mu, nu");

    // A inverse. A = V_perp + P_para
    TArray A;
    A("X, Y") = V_perp("X, Y") + P_para_("X, Y");
    auto A_eig = array_ops::array_to_eigen(A);
    using MatType = decltype(A_eig);
    MatType L_inv_eig = MatType(Eigen::LLT<MatType>(A_eig).matrixL()).inverse();
    auto tr0 = A.trange().data()[0];
    auto tr1 = A.trange().data()[1];
    assert(tr0 == tr1 && "Matrix A = LLT must be symmetric!");
    TArray L_inv =
        array_ops::eigen_to_array<Tile, Policy>(world, L_inv_eig, tr0, tr1);
    TArray A_inv;
    A_inv("X, Y") = L_inv("Z, X") * L_inv("Z, Y");

    // DF coeffs orthogonal to charge vector
    C_perp_("X, mu, nu") =
        A_inv("X, Y") * (Gamma_perp("Y, mu, nu") - Eta_perp("Y, mu, nu"));

    // total DF coeffs
    C_df_("X, mu, nu") = C_para_("X, mu, nu") + C_perp_("X, mu, nu");

    // build DF Coulomb term
    TArray J, J_part1, J_part2;
    J_part1("mu, nu") = 2.0 * Gamma_("X, mu, nu") * C_df_("X, rho, sig") * this->D_("rho, sig");
    J_part2("mu, nu") = C_df_("X, mu, nu") * V_("X, Y") * C_df_("Y, rho, sig") * this->D_("rho, sig");
    J("mu, nu") = J_part1("mu, nu") - J_part2("mu, nu");

    return J;
}

}  // namespace lcao

}  // namespace mpqc
