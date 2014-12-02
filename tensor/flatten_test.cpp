#include "../include/btas.h"
#include "../include/eigen.h"

#include <iostream>
#include <functional>

template <typename Extent>
std::vector<std::vector<unsigned long>>
data_order(Extent const &extent, unsigned long n) {
    auto ext_size = extent.size();

    std::vector<std::vector<unsigned long>> data_order;
    data_order.reserve(ext_size);

    // First dimension should always correspond to n
    std::vector<unsigned long> first_dim;
    first_dim.reserve(extent[n]);
    for (auto i = 0ul; i < extent[n]; ++i) {
        first_dim.push_back(i);
    }
    data_order.push_back(std::move(first_dim));

    // Remaining dims should go in order.
    for (auto dim = 0ul; dim < ext_size; ++dim) {
        if (dim != n) {
            std::vector<unsigned long> dim_vals;
            dim_vals.reserve(extent[dim]);
            for (auto elem = 0ul; elem < extent[dim]; ++elem) {
                dim_vals.push_back(elem);
            }
            data_order.push_back(std::move(dim_vals));
        }
    }

    return data_order;
}

template <std::size_t N, typename Extent>
Eigen::MatrixXd create_matrix(Extent const &extent) {
    auto rows = extent[N];
    auto cols = 1ul;
    for (auto i = 0ul; i < extent.size(); ++i) {
        if (i != N) {
            cols *= extent[i];
        }
    }
    return Eigen::MatrixXd(rows, cols);
}

// I would really like to use boost for this :(.
template <unsigned long Order, unsigned long N>
Eigen::MatrixXd matricization(btas::Tensor<double> const &input);

template <>
Eigen::MatrixXd matricization<3ul, 0ul>(btas::Tensor<double> const &input) {
    auto extent = input.extent();

    auto mat = create_matrix<0>(extent);

    for (auto i = 0ul; i < extent[0]; ++i) {
        for (auto j = 0ul; j < extent[1]; ++j) {
            for (auto k = 0ul; k < extent[2]; ++k) {
                mat(i, j *extent[2] + k) = input(i, j, k);
            }
        }
    }

    return mat;
}

template <>
Eigen::MatrixXd matricization<3ul, 1ul>(btas::Tensor<double> const &input) {
    auto extent = input.extent();

    auto mat = create_matrix<1>(extent);

    for (auto j = 0ul; j < extent[1]; ++j) {
        for (auto i = 0ul; i < extent[0]; ++i) {
            for (auto k = 0ul; k < extent[2]; ++k) {
                mat(j, i *extent[2] + k) = input(i, j, k);
            }
        }
    }

    return mat;
}

template <>
Eigen::MatrixXd matricization<3ul, 2ul>(btas::Tensor<double> const &input) {
    auto extent = input.extent();

    auto mat = create_matrix<2>(extent);

    for (auto k = 0ul; k < extent[2]; ++k) {
        for (auto i = 0ul; i < extent[0]; ++i) {
            for (auto j = 0ul; j < extent[1]; ++j) {
                mat(k, i *extent[1] + j) = input(i, j, k);
            }
        }
    }

    return mat;
}

template <std::size_t N>
void tucker_decomp(btas::Tensor<double> const &input); 

template <>
void tucker_decomp<3>(btas::Tensor<double> const &input) {
    auto rank = input.rank();

    std::vector<Eigen::MatrixXd> mats;
    mats.reserve(rank);

    mats.push_back(matricization<3,0>(input));
    mats.push_back(matricization<3,1>(input));
    mats.push_back(matricization<3,2>(input));

    for (auto const &mat : mats) {
        std::cout << "\nMat dims = " << mat.rows() << "x" << mat.cols() << "\n";
        std::cout << "Mat = \n" << mat << std::endl;
    }
    // TODO need to write KRON Mult and create low rank tensor type. 
}

int main() {
    btas::Tensor<double> tensor(2, 4, 3);

    auto i = 1;
    tensor.generate([=]() mutable { return i++; });

    Eigen::MatrixXd mat1(2, 12);
    Eigen::MatrixXd mat2(4, 6);
    Eigen::MatrixXd mat3(3, 8);


    for (auto i = 0; i < 2; ++i) {
        for (auto j = 0; j < 4; ++j) {
            for (auto k = 0; k < 3; ++k) {
                mat1(i, j * 3 + k) = tensor(i, j, k);
            }
        }
    }

    for (auto j = 0; j < 4; ++j) {
        for (auto i = 0; i < 2; ++i) {
            for (auto k = 0; k < 3; ++k) {
                mat2(j, i * 3 + k) = tensor(i, j, k);
            }
        }
    }

    for (auto k = 0; k < 3; ++k) {
        for (auto i = 0; i < 2; ++i) {
            for (auto j = 0; j < 4; ++j) {
                mat3(k, i * 4 + j) = tensor(i, j, k);
            }
        }
    }

    std::cout << "\nEigen Mat 1 = \n" << mat1 << std::endl;
    std::cout << "Eigen Mat 2 = \n" << mat2 << std::endl;
    std::cout << "Eigen Mat 3 = \n" << mat3 << std::endl;


    Eigen::JacobiSVD
        <Eigen::MatrixXd> svd(mat1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A1 = svd.matrixU().leftCols(2);

    svd.compute(mat2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A2 = svd.matrixU().leftCols(2);

    svd.compute(mat3, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd A3 = svd.matrixU().leftCols(2);

    Eigen::MatrixXd Kron(A3.cols() * A2.cols(), A3.rows() * A2.rows());

    for (auto i = 0; i < A3.cols(); ++i) {
        for (auto j = 0; j < A3.rows(); ++j) {
            for (auto k = 0; k < A2.cols(); ++k) {
                for (auto l = 0; l < A2.rows(); ++l) {
                    Kron(i * A2.cols() + k, j *A2.rows() + l) = A3(j, i)
                                                                * A2(l, k);
                }
            }
        }
    }


    Eigen::MatrixXd g_eig = A1.transpose() * mat1 * Kron.transpose();

    Eigen::MatrixXd mat1_approx = A1 * g_eig * Kron;
    std::cout << "Approximate mat1 = \n" << mat1_approx << std::endl;

    tucker_decomp<3>(tensor);

    return 0;
}
