#include <mpqc/chemistry/qc/scf/orbital_localization.h>

namespace mpqc {
namespace scf {

double boys_object(std::array<Mat, 3> const &xyz) {
    auto sum = 0.0;
    for (auto i = 0; i < xyz[0].cols(); ++i) {
        sum += xyz[0](i, i) * xyz[0](i, i);
        sum += xyz[1](i, i) * xyz[1](i, i);
        sum += xyz[2](i, i) * xyz[2](i, i);
    }

    return sum;
}

double gamma(double Aij, double Bij) {
    auto AB = std::sqrt(Aij * Aij + Bij * Bij);
    auto cos_gamma = -Aij / AB;
    auto sin_gamma = Bij / AB;
    auto ang = 0.25 * std::acos(cos_gamma) * ((sin_gamma < 0) ? -1 : 1);
    return (std::abs(ang) < 1e-7) ? 0 : ang;
};

void jacobi_sweeps(Mat &Cm, Mat &U, std::vector<Mat> const &ao_xyz) {
    std::array<Mat, 3> mo_xyz;
    mo_xyz[0] = Cm.transpose() * ao_xyz[0] * Cm;
    mo_xyz[1] = Cm.transpose() * ao_xyz[1] * Cm;
    mo_xyz[2] = Cm.transpose() * ao_xyz[2] * Cm;

    auto &mx = mo_xyz[0];
    auto &my = mo_xyz[1];
    auto &mz = mo_xyz[2];

    auto crit = boys_object(mo_xyz);
    auto iter = 1;
    auto error = crit - 0;
    while (error > 1e-4 && iter <= 50) {
        for (auto i = 0; i < Cm.cols(); ++i) {
            for (auto j = i + 1; j < Cm.cols(); ++j) {

                Vec3D vij = {mx(i, j), my(i, j), mz(i, j)};
                Vec3D vii = {mx(i, i), my(i, i), mz(i, i)};
                Vec3D vjj = {mx(j, j), my(j, j), mz(j, j)};

                double Aij = vij.squaredNorm()
                             - 0.25 * (vii - vjj).squaredNorm();
                double Bij = (vii - vjj).dot(vij);

                auto g = gamma(Aij, Bij);
                auto cg = std::cos(g);
                auto sg = std::sin(g);

                Eig::VectorXd col_Ui = U.col(i);
                Eig::VectorXd col_Uj = U.col(j);

                U.col(i) = cg * col_Ui + sg * col_Uj;
                U.col(j) = -sg * col_Ui + cg * col_Uj;

                for (auto z = 0; z < 3; ++z) {
                    auto &m = mo_xyz[z];
                    Eig::VectorXd z_i = m.col(i);
                    Eig::VectorXd z_j = m.col(j);

                    m.col(i) = cg * z_i + sg * z_j;
                    m.col(j) = -sg * z_i + cg * z_j;

                    z_i = m.row(i);
                    z_j = m.row(j);

                    m.row(i) = cg * z_i + sg * z_j;
                    m.row(j) = -sg * z_i + cg * z_j;
                }
            }
        }

        auto old_crit = crit;
        crit = boys_object(mo_xyz);
        error = std::abs(old_crit - crit);
        ++iter;
    }
}

} // namespace scf
} // namespace mpqc
