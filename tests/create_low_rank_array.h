#ifndef CREATE_LOW_RANK_ARRAY_H
#define CREATE_LOW_RANK_ARRAY_H

#include "../include/eigen.h"

namespace HeirChemTest {
template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
Matrix<T> low_rank_matrix(int rows, int cols, int rank) {
    auto full_rank = std::min(rows, cols);
    assert(rank <= std::min(rows, cols));

    Eigen::JacobiSVD<Matrix<T>> svd(Matrix<T>::Random(rows, cols),
                                    Eigen::ComputeThinU | Eigen::ComputeThinV);
    auto svals = svd.singularValues();

    for (auto i = 0;
         i < (full_rank - rank) && i < static_cast<int>(svals.size()); ++i) {
        svals[(svals.size() - 1) - i] = 0;
    }

    return svd.matrixU() * svals.asDiagonal() * svd.matrixV().transpose();
}

} // namespace HeirChemTest

#endif // CREATE_LOW_RANK_ARRAY_H
