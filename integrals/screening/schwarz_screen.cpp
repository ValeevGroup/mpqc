#include "schwarz_screen.h"

namespace mpqc {
namespace integrals {

Qmatrix::Qmatrix(MatrixD Q)
        : Q_(std::move(Q)),
          max_elem_in_row_(VectorD(Q_.rows())),
          max_elem_in_Q_(0.0)

{
    const auto nrows = Q_.rows();
    for (auto i = 0; i < nrows; ++i) {
        const auto max_elem = Q_.row(i).cwiseAbs().maxCoeff();
        max_elem_in_Q_ = std::max(max_elem, max_elem_in_Q_);
        max_elem_in_row_(i) = max_elem;
    }
}

} // namespace integrals
} // namespace mpqc
