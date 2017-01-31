#include "mpqc/chemistry/qc/lcao/integrals/screening/schwarz_screen.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

Qmatrix::Qmatrix(RowMatrixXd Q, std::unordered_map<int64_t, int64_t> map)
        : Q_(std::move(Q)),
          func_to_shell_map_(std::move(map)),
          max_elem_in_row_(Eigen::VectorXd(Q_.rows())),
          max_elem_in_Q_(0.0)

{
    const auto nrows = Q_.rows();
    for (auto i = 0; i < nrows; ++i) {
        const auto max_elem = Q_.row(i).cwiseAbs().maxCoeff();
        max_elem_in_Q_ = std::max(max_elem, max_elem_in_Q_);
        max_elem_in_row_(i) = max_elem;
    }
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
