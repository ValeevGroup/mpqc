#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_ALGEBRA_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_ALGEBRA_H_

#include "mpqc/math/external/eigen/eigen.h"

#include <madness/tensor/clapack.h>
#include <tiledarray.h>

#include "mpqc/math/tensor/clr/decomposed_tensor.h"

extern "C" void sgesdd_(const char *jobz, integer *m, integer *n, real4 *a,
                        integer *lda, real4 *s, real4 *u, integer *ldu,
                        real4 *vt, integer *ldvt, real4 *work, integer *lwork,
                        integer *iwork, integer *info);

extern "C" void dgesdd_(const char *jobz, integer *m, integer *n, real8 *a,
                        integer *lda, real8 *s, real8 *u, integer *ldu,
                        real8 *vt, integer *ldvt, real8 *work, integer *lwork,
                        integer *iwork, integer *info);

extern "C" void dgelqf_(integer *m, integer *n, real8 *a, integer *lda,
                        real8 *tau, real8 *work, integer *lwork,
                        integer *infoOUT);

extern "C" void dorglq_(integer *m, integer *n, integer *k, real8 *a,
                        integer *lda, real8 *tau, real8 *work, integer *lwork,
                        integer *info);

extern "C" void dpstrf_(const char *uplo, integer *n, real8 *a, integer *lda,
                        integer *piv, integer *rank, real8 *tol, real8 *work,
                        integer *info);

namespace mpqc {
namespace math {

static constexpr auto NoT = madness::cblas::CBLAS_TRANSPOSE::NoTrans;
static constexpr auto Tr = madness::cblas::CBLAS_TRANSPOSE::Trans;

// Compute the column pivoted qr decomposition into data, will modify input
// pointers data and J
integer col_pivoted_qr(double *data, double *Tau, integer rows, integer cols,
                       integer *J);

integer non_pivoted_qr(double *data, double *Tau, integer rows, integer cols);

integer non_pivoted_lq(double *data, double *Tau, integer rows, integer cols);

integer svd(double *data, double *s, double *u, double *vt, integer rows,
            integer cols);

integer svd(double *data, double *s, double *u, double *vt, integer rows,
            integer cols, const char JOBZ);

integer form_q(double *data, double *Tau, integer rows, integer rank);

integer form_q_from_lq(double *data, double *Tau, integer cols, integer rows,
                       integer rank);

size_t svd_rank(double const *s, size_t N, double thresh);

// calculate the qr_rank of a matrix.
std::size_t qr_rank(double const *data, std::size_t rows, std::size_t cols,
                    double threshold);

// Returns true if input is low rank.
bool full_rank_decompose(
    TA::Tensor<double> const &in, TA::Tensor<double> &L, TA::Tensor<double> &R,
    double thresh,
    std::size_t max_out_rank = std::numeric_limits<std::size_t>::max());

void ta_tensor_col_pivoted_qr(TA::Tensor<double> &in, TA::Tensor<double> &L,
                              TA::Tensor<double> &R, double thresh);

// eats in data and outputs L and R tensors.
void ta_tensor_qr(TA::Tensor<double> &in, TA::Tensor<double> &L,
                  TA::Tensor<double> &R);

// eats in data and outputs L and R tensors.
void ta_tensor_lq(TA::Tensor<double> &in, TA::Tensor<double> &L,
                  TA::Tensor<double> &R);

// Sacrifice input data
void ta_tensor_svd(TA::Tensor<double> &in, TA::Tensor<double> &L,
                   TA::Tensor<double> &R, double thresh);

/// Currently modifies input data regardless could cause some loss of accuracy.
void recompress(DecomposedTensor<double> &t);

/// Returns an empty DecomposedTensor if the compression rank was to large.
DecomposedTensor<double> two_way_decomposition(
    DecomposedTensor<double> const &t);

TA::Tensor<double> combine(DecomposedTensor<double> const &t);

/*! \brief performs the pivoted Cholesky decomposition
 *
 * This function will take a matrix that is symmetric semi-positive definate
 * and over write it with the permuted and truncated Cholesky vectors
 * corresponding to the matrix.
*/
integer piv_cholesky(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> &a);

}  // namespace math
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_ALGEBRA_H_
