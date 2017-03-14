//
// Created by Chong Peng on 2/22/17.
//

#ifndef SRC_MPQC_MATH_LINALG_GRAM_SCHMIDT_H_
#define SRC_MPQC_MATH_LINALG_GRAM_SCHMIDT_H_

#include <vector>

#include <TiledArray/algebra/utils.h>

/**
 * Gram-Schmidt Process
 *
 * reference:
 * https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Algorithm
 */

template <typename D>
void gram_schmidt(std::vector<D>& V, std::size_t start = 0) {
  const auto k = V.size();

  TA_ASSERT(start < k);

  for (std::size_t i = start; i < k; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      typename D::element_type tmp = dot_product(V[i], V[j]);
      axpy(V[i], -tmp, V[j]);
    }
    // normalize
    scale(V[i], 1.0 / norm2(V[i]));
  }

  // test
  //#ifndef NDEBUG
  //  const auto tolerance = std::numeric_limits<typename
  //  D::element_type>::epsilon()*100;
  //  for(auto i = 0; i < k; ++i){
  //    for(auto j = i ; j < k; ++j ){
  //      const auto test = dot_product(V[i], V[j]);
  //      std::cout << "i= " << i << " j= " << j << " dot= " << test <<
  //      std::endl;
  //      if(i==j){
  //        TA_ASSERT( test - 1.0 < tolerance);
  //      }
  //      else{
  //        TA_ASSERT( test < tolerance);
  //      }
  //    }
  //  }
  //#endif
}


/**
 *  vector V1 is already orthonormalized
 *  orthonormalize V2 with respect to V1
 */
template <typename D>
void gram_schmidt(const std::vector<D>& V1, std::vector<D>& V2) {
  const auto k1 = V1.size();
  const auto k2 = V2.size();

  for (std::size_t i = 0; i < k2; ++i) {
    // loop over all vector in V1
    for (std::size_t j = 0; j < k1; ++j) {
      typename D::element_type tmp = dot_product(V2[i], V1[j]);
      axpy(V2[i], -tmp, V1[j]);
    }

    // loop over other vector in V2
    for (std::size_t k = 0; k < i; ++k) {
      typename D::element_type tmp = dot_product(V2[i], V2[k]);
      axpy(V2[i], -tmp, V2[k]);
    }

    // normalize
    scale(V2[i], 1.0 / norm2(V2[i]));
  }
}

#endif  // SRC_MPQC_MATH_LINALG_GRAM_SCHMIDT_H_
