#ifndef MPQC4_SRC_MPQC_MATH_CP_ALS_H
#define MPQC4_SRC_MPQC_MATH_CP_ALS_H

#ifdef TILEDARRAY_HAS_BTAS

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <TiledArray/config.h>

#include <TiledArray/conversions/btas.h>
#include <TiledArray/external/btas.h>
#include <TiledArray/tiled_range.h>
#ifdef HAVE_INTEL_MKL
#define _HAS_INTEL_MKL
#endif
#include <btas/generic/cp_als.h>
#include <btas/generic/cp_rals.h>
namespace mpqc {
namespace math {

/**
 * REQUIRES BTAS
 * Calculates the canonical product (CP) decomposition of /c reference using
 * alternating least squares (ALS) method in BTAS
 *
 * @param[in,out] reference In: tensor to be decomposed. Out: if \c recompose
 * the CP approximated tensor else unchanged
 * @param[in, out] factor_matrices In: empty vector. Out: Set of N CP-factor
 * matrices.
 * @param[in] recompose false Return the CP approximation of tensor \c reference
 * @param[in] optimized_rank true Compute the CP decomposition to finite error
 * (i.e. best rank decomposition of \c reference)
 * @param[in] geometric false If \c optimized_rank = false build \c rank
 * approximation using geometric step size
 * @param[in] tcutCP 0.01 if \c optimized_rank = true, error threshold for
 * difference in \c reference and its approximation
 * @param[in] rank 0 if \c optimized_rank = false, what is the CP decomposition
 * rank
 * @param[in] direct true Should the method compute the ALS without the
 * Khatri-Rao product
 * @param[in] calculate_epsilon fasle if \c optimized_rank = false should the
 * method compute the error between CP approximation and \c reference
 * @param[in] geometric_step 1 if \c geometric = true what is the step size
 * @param[in] step 1 approximation is built from r=1 to r = (\c optimized_rank)
 * ? R : \c rank with step size of \c step
 * @param[in] max_rank 1e4 if \c optimized_rank = true max computed rank, R
 * @param[in] max_als 1e3 max number of iterations to optimize rank r factors
 * using ALS
 * @param[in] tcutALS 0.1 ALS error threshold
 * @param[in] SVD_initial_guess false Should the initial guess be computed using
 * the from the left singular vectors
 * @param[in] SVD_rank 0 if \c SVD_initial_guess = true how many left singular
 * vectors (rank of initial guess) \c SVD_rank <= smallest dimension of \c
 * reference.
 */

;

template <typename Array>
void cp_rals_compute_rank(Array &reference, std::vector<Array> & factor_matrices, int block_size = 1, bool recompose = false,
                          int rank = 1, bool direct = true, bool calculate_epsilon = false, int step = 1, int max_als = 1000,
                          double ALSThresh = 0.1, bool SVD_initial_guess = false, int SVD_rank = 0, bool symm = false);

    template <typename Array>
    void cp_als(Array &reference, std::vector<Array> &factor_matrices, int block_size = 1,
                bool recompose = false, bool optimized_rank = true,
                bool geometric = false, double tcutCP = 0.01, int rank = 0,
                bool direct = true, bool calculate_epsilon = false,
                int geometric_step = 1, int step = 1, int max_rank = 1e4,
                double max_als = 1e3, double tcutALS = 0.1,
                bool SVD_initial_guess = false, int SVD_rank = 0, bool symm = false) {
  madness::World &world = reference.world();
  auto one_node = (world.size() == 1);

  // get the array information and transform it to a btas object
  using Tile = typename Array::value_type;
  using Policy = typename Array::policy_type;

  // transform the array to a btas object
  reference.make_replicated();
  world.gop.fence();

  auto btas_array = TiledArray::array_to_btas_tensor<Tile, Policy>(reference);

  std::vector<decltype(btas_array)> btas_factors_matrices;
  int ndim = btas_array.rank();

  if (world.rank() == 0) {
    // compute the CP decomposition based on user's interest
    btas::CP_ALS<decltype(btas_array)> CP(btas_array);
    auto error = 0.0;
    if (optimized_rank) {
      CP.compute_error(tcutCP, direct, step, max_rank, max_als, tcutALS,
                               SVD_initial_guess, SVD_rank);
    }
    else if (geometric) {
      CP.compute_geometric(rank, geometric_step, direct, max_als,
                           calculate_epsilon, tcutALS, SVD_initial_guess,
                           SVD_rank);
    }
    else {
      error = CP.compute_rank(rank, direct, calculate_epsilon, step, max_als, tcutALS,
                      SVD_initial_guess, SVD_rank, symm);
    }

    //std::cout << "Error in CP decomposition is " << error << "\nnorm is " <<  TiledArray::norm2(reference) << std::endl;
    // obtain factor matrices prepare to push them back as TA objects
    btas_factors_matrices = CP.get_factor_matrices();
    int brank = btas_factors_matrices[0].extent(1);

    if (recompose) {
      btas_array = CP.reconstruct();
    }

    // Scale the first factor matrix by the parallel factor, this choice is
    // arbitrary
    for (int i = 0; i < brank; i++) {
      btas::scal(btas_factors_matrices[0].extent(0),
                 btas_factors_matrices[ndim](i),
                 std::begin(btas_factors_matrices[0]) + i, brank);
    }
  }
  
  // If one wants to reconstruct the full tensor from the factor matrices it
  // will be stored in the reference passed in.
  if (recompose) {
    // Share btas_array to all nodes
    world.gop.broadcast_serializable(btas_array, 0);

    // Convert btas_array back to a TiledArray object
    auto trange = reference.trange();
    reference = TiledArray::btas_tensor_to_array<Array>(world, trange,
                                                        btas_array, !one_node);
  }

  // Take each btas factor matrix turn it into a TA factor matrix.
  // Fill the vector which stores the TA factor matrices
  // Share the BTAS tensor vector of factor matrices
  world.gop.broadcast_serializable(btas_factors_matrices, 0);
  for (auto factor = 0; factor < ndim; factor++) {
    auto row_trange = reference.trange().data()[factor];
    auto col_trange = mpqc::utility::compute_trange1(btas_factors_matrices[0].extent(1), block_size);
    //TiledArray::TiledRange1 col_trange(0, btas_factors_matrices[0].extent(1));
    TiledArray::TiledRange trange({row_trange, col_trange});

    auto TA_factor = TiledArray::btas_tensor_to_array<Array>(
        world, trange, btas_factors_matrices[factor], !one_node);
    factor_matrices.push_back(TA_factor);
  }

  return;
}

template <typename Array>
void cp_rals_compute_rank(Array &reference, std::vector<Array> & factor_matrices, int block_size, bool recompose,
                          int rank, bool direct, bool calculate_epsilon, int step, int max_als,
                          double ALSThresh, bool SVD_initial_guess, int SVD_rank, bool symm){
  auto & world = reference.world();
  auto one_node = (world.size() == 1);

  using Tile = typename Array::value_type;
  using Policy = typename Array::policy_type;

  // transform the array to a btas object
  reference.make_replicated();

  auto btas_array = TiledArray::array_to_btas_tensor<Tile, Policy>(reference);
  std::vector<decltype(btas_array)> btas_factors_matrices;
  int ndim = btas_array.rank();

  if(world.rank() == 0){
    //auto t1 = std::chrono::high_resolution_clock::now();
    btas::CP_RALS<decltype(btas_array)> CP(btas_array);
    auto diff = CP.compute_rank(rank, direct, calculate_epsilon, step, max_als, ALSThresh,
                      SVD_initial_guess, SVD_rank, symm);
    //auto t2 = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double> time_span = t2 - t1;
    //std::cout << "Time to compute cp-rals is " << time_span.count() << " diff is " << diff << std::endl;
    // obtain factor matrices prepare to push them back as TA objects
    btas_factors_matrices = CP.get_factor_matrices();
    int brank = btas_factors_matrices[0].extent(1);

    if (recompose) {
      btas_array = CP.reconstruct();
    }

    // Scale the first factor matrix by the parallel factor, this choice is
    // arbitrary
    for (int i = 0; i < brank; i++) {
      btas::scal(btas_factors_matrices[0].extent(0),
                 btas_factors_matrices[ndim](i),
                 std::begin(btas_factors_matrices[0]) + i, brank);
    }
  }

  // If one wants to reconstruct the full tensor from the factor matrices it
  // will be stored in the reference passed in.
  if (recompose) {
    // Share btas_array to all nodes
    world.gop.broadcast_serializable(btas_array, 0);

    // Convert btas_array back to a TiledArray object
    auto trange = reference.trange();
    reference = TiledArray::btas_tensor_to_array<Array>(world, trange,
                                                        btas_array, !one_node);
  }

  world.gop.broadcast_serializable(btas_factors_matrices, 0);
  for (auto factor = 0; factor < ndim; factor++) {
    auto row_trange = reference.trange().data()[factor];
    //TiledArray::TiledRange1 col_trange(0, btas_factors_matrices[0].extent(1));
    auto col_trange = mpqc::utility::compute_trange1(btas_factors_matrices[0].extent(1), block_size);
    TiledArray::TiledRange trange({row_trange, col_trange});

    auto TA_factor = TiledArray::btas_tensor_to_array<Array>(
        world, trange, btas_factors_matrices[factor], !one_node);
    factor_matrices.push_back(TA_factor);
  }
  return;
}

}  // namespace math
}  // namespace mpqc
#endif  // TILEDARRAY_HAS_BTAS
#endif  // MPQC4_SRC_MPQC_MATH_CP_ALS_H
