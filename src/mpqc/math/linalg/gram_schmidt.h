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
void gram_schmidt(std::vector<D>& V, std::size_t start=0) {
  const auto k = V.size();

  TA_ASSERT(start < k);

  for (auto i = start; i < k; ++i) {
    for (auto j = 0; j < i; ++j) {
      typename D::element_type tmp = dot_product(V[i], V[j]);
      axpy(V[i], -tmp, V[j]);
    }
    // normalize
    scale(V[i], 1.0 / norm2(V[i]));
  }

  // test
//#ifndef NDEBUG
//  const auto tolerance = std::numeric_limits<typename D::element_type>::epsilon()*100;
//  for(auto i = 0; i < k; ++i){
//    for(auto j = i ; j < k; ++j ){
//      const auto test = dot_product(V[i], V[j]);
//      std::cout << "i= " << i << " j= " << j << " dot= " << test << std::endl;
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

#endif  // SRC_MPQC_MATH_LINALG_GRAM_SCHMIDT_H_
