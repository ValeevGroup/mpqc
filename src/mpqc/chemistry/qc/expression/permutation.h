//
// Created by Chong Peng on 8/16/16.
//

#ifndef MPQC_PERMUTATION_H
#define MPQC_PERMUTATION_H

#include "formula.h"

namespace mpqc{

namespace detail{

  /// swap bra with ket
  void swap(Formula& formula);

  /// swap inside bra or ket
  void swap_internal(Formula& formula, bool right_size=true);

  /// swap between bra and ket
  void swap_external(Formula& formula, bool right_size=true);

  /// return unique permutation of ( p q | r s )
  std::vector<Formula> permutations_chemical(const Formula& formula);

  /// return unique permutation of < p q | r s >
  std::vector<Formula> permutations_physical(const Formula& formula);

}

std::vector<Formula> permutations(const Formula& formula);

}// end of namespace mpqc


#endif //MPQC_PERMUTATION_H
