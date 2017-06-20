//
// Created by Chong Peng on 10/15/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_FORMULA_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_FORMULA_H_

#include <vector>

#include "mpqc/math/groups/petite_list.h"

#include "greek_to_english_name.h"
#include "mpqc/chemistry/qc/lcao/expression/operator.h"
#include "orbital_index.h"

#include <iostream>

using mpqc::lcao::OrbitalIndex;
using mpqc::Operator;

namespace mpqc {

namespace detail {
/// provides identity transform for strings
struct identity {
  std::string operator()(std::string &&arg) const {
    return std::forward<std::string>(arg);
  }
};
/// transforms strings by appending string representation of an integral counter
struct append_count {
  append_count(long count_init = 0) : count_(count_init) {}
  append_count(const append_count &append_count) = default;
  std::string operator()(std::string &&arg) {
    return arg + std::to_string(count_++);
  }
  long count_;
};
}

/**
 * \brief Formula parses a string represnetation of quantum mechanical matrix elements and
 *        related expressions.
 */
class Formula {
 public:
  /// Types of Notation
  enum class Notation { Invalid = -1, Chemical = 0, Physical = 1 };
  /// Position = Bra or Ket
  enum class Position { Invalid = -1, Bra = 0, Ket = 1 };

  /// Types of Options
  enum class Option { DensityFitting = 0, Inverse = 1, InverseSquareRoot = 2 };

  static const std::map<Option, std::wstring> option_to_string;

  Formula() {}
  Formula(Formula const &) = default;
  Formula(Formula &&) = default;
  Formula &operator=(Formula const &) = default;
  Formula &operator=(Formula &&) = default;

  /**
 *  Constructor parses the string in one of the following formats:
 *  - Physical Notation
 *  <index1 index2|operator|index3 index4>[option1,option2]
 *
 *  - Chemical Notation
 *  (index1 index2|operator|index3 index4)[option]
 *
 *  where the index keys \c index1 , \c \index2 , etc. are parsed by ::mpqc::lcao::OrbitalIndex,
 *  the operator key \c operator
 *  is parsed by Operator , and the option keys can be one of the following:
 *  - \c df -> Option::DensityFitting
 *  - \c inv -> Option::Inverse
 *  - \c inv_sqr -> Option::InverseSquareRoot
 *  - \c <symm> -> math::PetiteList::Symmetry
 */
  Formula(std::wstring formula);

  /// reconstruct a std::wstring representation of the formula
  /// @sa Formula::to_ta_expression()
  std::wstring string() const;

  /// @name index functions

  /// dimension of formula(2, 3 or 4)
  std::size_t rank() const;

  /// @return true if this contains \c index
  bool has_index(const OrbitalIndex &index) const;

  /// @return true if it only contains AO indices
  bool is_ao() const;

  /// return bra_index
  const std::vector<OrbitalIndex> &bra_indices() const { return bra_indices_; }

  /// return bra_index
  std::vector<OrbitalIndex> &bra_indices() { return bra_indices_; }

  /// return ket_index
  const std::vector<OrbitalIndex> &ket_indices() const { return ket_indices_; }

  /// return ket_index
  std::vector<OrbitalIndex> &ket_indices() { return ket_indices_; }

  /// set ket_index
  void set_ket_indices(const std::vector<OrbitalIndex> &bra_idxs) {
    ket_indices_ = bra_idxs;
  }

  /// set bra index
  void set_bra_indices(const std::vector<OrbitalIndex> &ket_idxs) {
    bra_indices_ = ket_idxs;
  }
  /// @}

  /// @name Notation functions
  /// @{

  /// Notation accessor
  const Notation &notation() const {
    TA_USER_ASSERT(notation_ != Notation::Invalid, "invalid Notation")
    return notation_;
  }

  /// set Notation
  void set_notation(const Notation &notation) {
    TA_USER_ASSERT(notation != Notation::Invalid, "invalid Notation")
    notation_ = notation;
  }

  /// @}

  /// @name Operator functions
  /// @{

  /// set Operator
  void set_operator(const Operator &oper) { oper_ = oper; }

  /// set Operator type @sa Operator::Type
  void set_operator_type(const Operator::Type &oper_type) {
    oper_.set_type(oper_type);
  }

  /// Operator accessor
  const Operator &oper() const { return oper_; }

  /// @}

  /// Symmetry accessor
  math::PetiteList::Symmetry symmetry() const {
    return symm_;
  }

  /// @name Formula options functions
  /// @{

  /// option accessor
  const std::vector<Option>& options() const {
    return options_;
  }

  /// @param op a Formula::Option object
  /// Calling this ensures that \c has_option(op) will return \c true
  void add_option(Option op);

  /// @param op a Formula::Option object
  /// @return true if this formula has option \c op
  bool has_option(Option op) const;

  /// clear vector<Option>
  void clear_option();

  /// remove Option op in options_
  void remove_option(Option op);

  /// @}

  /// @name Comparison operators
  /// @{

  bool operator<(const Formula &other) const;
  bool operator==(const Formula &other) const;
  bool operator!=(const Formula &other) const { return !(*this == other); }

  /// @}

  /// converts this to a TA expression annotation
  /// @tparam Transformer a unary functor class
  /// @param transform_op used to transform index keys
  template <typename Transformer = detail::identity>
  std::string to_ta_expression(
      Transformer transform_op = detail::identity()) const;

 private:
  /// parse the index on one side
  std::vector<OrbitalIndex> check_orbital_index(std::wstring index_array);

 private:
  Operator oper_;
  Notation notation_ = Notation::Invalid;
  std::vector<Option> options_;
  math::PetiteList::Symmetry symm_ = math::PetiteList::Symmetry::e;
  std::vector<OrbitalIndex> bra_indices_;
  std::vector<OrbitalIndex> ket_indices_;
};


template <typename Transformer>
std::string Formula::to_ta_expression(Transformer transform_op) const {
  std::string ta_expression;
  std::size_t rank = this->rank();
  std::size_t count = 0;

  // add left index
  for (const auto &index : bra_indices_) {
    std::string index_expression = transform_op(index.to_ta_expression());
    ta_expression.append(index_expression.begin(), index_expression.end());
    ++count;
    ta_expression.append(", ");
  }

  // add right index
  for (const auto &index : ket_indices_) {
    std::string index_expression = transform_op(index.to_ta_expression());
    ta_expression.append(index_expression.begin(), index_expression.end());
    ++count;
    if (count < rank) {
      ta_expression.append(", ");
    }
  }

  return ta_expression;
}
} // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_FORMULA_H_
