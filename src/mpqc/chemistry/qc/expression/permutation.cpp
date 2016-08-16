//
// Created by Chong Peng on 8/16/16.
//

#include "permutation.h"


namespace mpqc{

namespace detail{

void swap(Formula& formula)
{

  auto bra = formula.bra_indices();
  formula.set_bra_indices(formula.ket_indices());
  formula.set_ket_indices(bra);

}

void swap_internal(Formula& formula, bool right_size)
{
  if(right_size){
    auto ket = formula.ket_indices();
    std::swap(ket[0],ket[1]);
    formula.set_ket_indices(ket);
  }else{
    auto bra = formula.bra_indices();
    std::swap(bra[0],bra[1]);
    formula.set_bra_indices(bra);
  }
}

void swap_external(Formula& formula, bool right_size){
  // swap index 1
  if(right_size){
    auto index = formula.bra_indices()[1];
    formula.bra_indices()[1] = formula.ket_indices()[1];
    formula.ket_indices()[1] = index;
  }
  // swap index 0
  else{
    auto index = formula.bra_indices()[0];
    formula.bra_indices()[0] = formula.ket_indices()[0];
    formula.ket_indices()[0] = index;
  }

}

std::vector<Formula> permutations_chemical(const Formula& formula){
  std::vector<Formula> result;
  auto bra = formula.bra_indices();
  auto ket = formula.ket_indices();

  if(bra[0]!=bra[1]){
    Formula permutation = formula;
    swap_internal(permutation,false);
    result.push_back(permutation);
  }

  if(ket[0]!=ket[1]){
    Formula permutation = formula;
    swap_internal(permutation,true);
    result.push_back(permutation);
  }

  if( (bra[0]!=bra[1]) && (ket[0]!=ket[1]) ){
    Formula permutation = formula;
    swap_internal(permutation,true);
    swap_internal(permutation,false);
    result.push_back(permutation);
  }

  return result;
}

std::vector<Formula> permutations_physical(const Formula& formula){
  std::vector<Formula> result;
  auto bra = formula.bra_indices();
  auto ket = formula.ket_indices();

  if(bra[0]!=ket[0]){
    Formula permutation = formula;
    swap_external(permutation,false);
    result.push_back(permutation);
  }

  if(bra[1]!=ket[1]){
    Formula permutation = formula;
    swap_external(permutation,true);
    result.push_back(permutation);
  }

  if( (bra[0]!=ket[0]) && (bra[1]!=ket[1]) ){
    Formula permutation = formula;
    swap_external(permutation,true);
    swap_external(permutation,false);
    result.push_back(permutation);
  }

  return result;
}

} // end of namespace detail

std::vector<Formula> permutations(const Formula& formula){
  std::vector<Formula> result;

  if(formula.rank()==2){
    Formula permutation = formula;
    detail::swap(permutation);
    result.push_back(permutation);
  }
  else if(formula.rank()==3){

    Formula permutation = formula;
    detail::swap_internal(permutation);
    result.push_back(permutation);
  }
  else if(formula.rank()==4){

    //chemical notation
    if(formula.notation() == Formula::Notation::Chemical){

      auto permutations1 = detail::permutations_chemical(formula);

      result.insert(result.cend(),permutations1.begin(), permutations1.end());

      if(formula.bra_indices() != formula.ket_indices()){
        Formula swap_bra_ket = formula;
        detail::swap(swap_bra_ket);

        auto permutations2 = detail::permutations_chemical(formula);
        result.insert(result.cend(),permutations2.begin(), permutations2.end());
      }

    }
    //physical notation
    else if(formula.notation() == Formula::Notation::Physical){

      auto permutations1 = detail::permutations_physical(formula);

      result.insert(result.cend(),permutations1.begin(), permutations1.end());

      if(formula.bra_indices() != formula.ket_indices()){
        Formula swap_bra_ket = formula;
        detail::swap(swap_bra_ket);

        auto permutations2 = detail::permutations_physical(formula);
        result.insert(result.cend(),permutations2.begin(), permutations2.end());
      }

    }
    else{
      throw  std::runtime_error("Invalid Notation!");
    }

  }

  return result;

}

} // end of namespace mpqc
