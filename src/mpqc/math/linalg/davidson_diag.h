//
// Created by ChongPeng on 2/20/17.
//

#ifndef MPQC_DAVIDSON_DIAG_H
#define MPQC_DAVIDSON_DIAG_H

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/util/misc/exenv.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"

#include <tiledarray.h>
#include <TiledArray/algebra/utils.h>
#include <mpqc/util/misc/assert.h>

namespace mpqc {

namespace detail {

// template <typename Tile, typename Policy>
// TA::DistArray<Tile, Policy> matrix_multiply(
//    const TA::DistArray<Tile, Policy>& left,
//    const TA::DistArray<Tile, Policy>& right) {
//  TA::DistArray<Tile, Policy> result;
//  const auto left_rank = left.size();
//  const auto right_rank = right.size();
//
//  if (left_rank == 2 && right_rank == 2) {
//    result("i,j") = left("i,k") * right("k,j");
//  } else {
//    ExEnv::out0() << "left rank = " << left_rank
//                  << " right rank = " << right_rank << "\n";
//    throw ProgrammingError("matrix_multiply with this rank not
//    implemented!\n",
//                           __FILE__, __LINE__);
//  }
//  return result;
//};

template <typename Tile, typename Policy>
RowMatrix<typename TA::DistArray<Tile, Policy>::value_type> dot_product(
    const TA::DistArray<Tile, Policy>& left,
    const TA::DistArray<Tile, Policy>& right) {
  const auto left_rank = left.size();
  const auto right_rank = right.size();
  TA::DistArray<Tile, Policy> result;

  if (left_rank == 2 && right_rank == 2) {
    result("i,j") = left("i,k") * right("k,j");
  } else {
    ExEnv::out0() << "left rank = " << left_rank
                  << " right rank = " << right_rank << "\n";
    throw ProgrammingError("dot_product with this rank not implemented !\n ",
                           __FILE__, __LINE__);
  }

  return array_ops::array_to_eigen(result);
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> transpose(
    const TA::DistArray<Tile, Policy>& array) {
  TA::DistArray<Tile, Policy> result;
  const auto rank = array.size();

  if (rank == 2) {
    result("i,j") = array("j,i");

  } else {
    ExEnv::out0() << "rank = " << rank << "\n";
    throw ProgrammingError("transpose with this rank not implemented!\n",
                           __FILE__, __LINE__);
  }
  return result;
};

template <typename Tile, typename Policy>
void eigen_to_datatype(
    const RowMatrix<typename TA::DistArray<Tile, Policy>::value_type>&
        eigen_matrix,
    TA::DistArray<Tile, Policy>& ta_array) {
  const auto cols = eigen_matrix.cols();
  const auto rows = eigen_matrix.rows();
  // block this by one tile
  TA::TiledRange1 col_tr = {0, cols};
  TA::TiledRange1 row_tr = {0, rows};

  ta_array = array_ops::eigen_to_array<Tile, Policy>(
      TA::get_default_world(), eigen_matrix, col_tr, row_tr);
};

/*
 * use QR to orthognolize two rank 2 TA::DistArray
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> orthognolize(
    const TA::DistArray<Tile, Policy>& left,
    const TA::DistArray<Tile, Policy>& right) {
  const auto left_rank = left.size();
  const auto right_rank = right.size();

  const auto rows_trange1 = left.trange()[0];
  // assert same number of cols
  TA_ASSERT(rows_trange1.extent() == right.trange()[0].extent());

  if (left_rank == 2 && right_rank == 2) {
    auto left_trange1 = left.trange()[1];
    auto right_trange1 = right.trange()[1];

    auto left_eigen = array_ops::array_to_eigen(left);
    auto right_eigen = array_ops::array_to_eigen(right);

    RowMatrix<typename TA::DistArray<Tile, Policy>::element_type>
        left_merge_right(left_eigen.rows(),
                         left_eigen.cols() + right_eigen.cols());
    left_merge_right << left_eigen, right_eigen;

    // perform QR
    Eigen::HouseholderQR<decltype(left_merge_right)> qr(left_merge_right);

    decltype(left_merge_right) Q = qr.householderQ();

    // convert Q into TA::DistArray
    auto left_merge_right_trange1 =
        utility::detail::join_trange1(left_trange1, right_trange1);

    auto result =
        array_ops::eigen_to_array(TA::get_default_world(), left_merge_right,
                                  rows_trange1, left_merge_right_trange1);

    return result;

  } else {
    ExEnv::out0() << "left rank = " << left_rank
                  << " right rank = " << right_rank << "\n";
    throw ProgrammingError("orthognolize with this rank not implemented !\n ",
                           __FILE__, __LINE__);
  }

};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> davidson_diag_residual(
    TA::DistArray<Tile, Policy>& HB, TA::DistArray<Tile, Policy>& B,
    TA::DistArray<Tile, Policy>& V,
    const RowVector<typename TA::DistArray<Tile, Policy>::element_type>& E){

  const auto HB_rank = HB.size();
  const auto B_rank = B.size();

  if(HB_rank == 2 and B_rank == 2){
      TA::DistArray<Tile,Policy> BV;
      BV("i,j") = B("i,k")*V("k,j");
      HB("i,j") = HB("i,k")*BV("k,j");

      // lambda to update BV
      auto update = [&E] (Tile& result_tile) {

        const auto& range = result_tile.range();
        double norm = 0.0;
        for (const auto& i : range){
          const auto result = result_tile[i] * E[i[1]];
          result_tile[i] = result;
          norm += result*result;
        }
        return std::sqrt(norm);
      };

      TA::foreach_inplace(BV, update);
      B.world().gop.fence();
      HB("i,j") -= BV("i,j");

      return HB;

  }else{
    ExEnv::out0() << "HB rank = " << HB_rank
                  << " B rank = " << B_rank << "\n";
    throw ProgrammingError("davidson_diag_residual with this rank not implemented !\n ",
                           __FILE__, __LINE__);
  }

};

}  // namespace detail

/**
 * \brief Davidson Algorithm
 *
 * \tparam D array type
 *
 *
 */

template <typename D>
class SymmDavidsonDiag {
 public:
  using value_type = typename D::element_type;
  using result_type = RowVector<value_type>;

  SymmDavidsonDiag(unsigned int n_roots, unsigned int n_guess)
      : n_roots_(n_roots), n_guess_(n_guess) {}

  RowVector<value_type> extrapolate(D& HB, D& B) {
    // subspace
    // dot_product will return a replicated Eigen Matrix
    RowMatrix<value_type> G = detail::dot_product(detail::transpose(B), HB);

    // do eigen solve locally
    result_type E;
    D V;
    {
      // this return eigenvalue and eigenvector with size n_guess
      Eigen::EigenSolver<RowMatrix<value_type>> es(G);
      RowMatrix<value_type> v = es.eigenvectors().real().leftCols(n_guess_);
      E = es.real().eigenvalues().segment(0, n_guess_);
      detail::eigen_to_datatype(v, V);
    }

    // compute residual
    D B_new = detail::davidson_diag_residual(HB, B, V, E);

    // orthognolize B_new and update B
    B = detail::orthognolize(B, B_new);

    return E.segment(0, n_roots_);
  }

 private:
  unsigned int n_roots_;
  unsigned int n_guess_;
};

}  // namespace mpqc

#endif  // MPQC_DAVIDSON_DIAG_H
