//
// Created by Chong Peng on 10/15/15.
//

#ifndef MPQC_FORMULA_H
#define MPQC_FORMULA_H

#include <vector>

#include "greek_to_english_name.h"
#include "operator.h"
#include "orbital_index.h"

using mpqc::OrbitalIndex;
using mpqc::Operator;

namespace mpqc {

namespace detail {
/// provides identity transform for strings
struct identity {
  std::string operator()(std::string&& arg) const {
    return std::forward<std::string>(arg);
  }
};
/// transforms strings by appending string representation of an integral counter
struct append_count {
  append_count(long count_init = 0) : count_(count_init) {}
  append_count(const append_count& append_count) = default;
  std::string operator()(std::string&& arg) {
    return arg + std::to_string(count_++);
  }
  long count_;
};
}

/**
 * \brief Formula class that represent quantum mechanical expressions
 *
 * format for formula
 *
 *  - Physical Notation
 *  <index1 index2|operation|index3 index4>[option1,option2]
 *
 *  - Chemical Notation
 *  (index1 index2|operation|index3 index4)[option]
 *
 *  @sa OrbitalIndex for description of index
 *  @sa Operation for description of operation and option
 */
class Formula {
 public:
  /// Types of Notation
  enum class Notation { Invalid = -1, Chemical = 0, Physical = 1 };
  /// Position = Bra or Ket
  enum class Position { Invalid = -1, Bra = 0, Ket = 1 };

  Formula() : notation_(Notation::Invalid) {}
  Formula(Formula const &) = default;
  Formula(Formula &&) = default;
  Formula &operator=(Formula const &) = default;
  Formula &operator=(Formula &&) = default;

  /**
   *  Constructor
   *  @param formula a properly formatted std::wstring
   */
  Formula(std::wstring formula);

  /// reconstruct a std::wstring representation of the formula
  /// @sa Formula::to_ta_expression()
  std::wstring string() const;

  /// return left_index
  const std::vector<OrbitalIndex> &bra_indices() const { return bra_indices_; }

  /// return left_index
  std::vector<OrbitalIndex> &bra_indices() { return bra_indices_; }

  /// return right_index
  const std::vector<OrbitalIndex> &ket_indices() const { return ket_indices_; }

  /// return right_index
  std::vector<OrbitalIndex> &ket_indices() { return ket_indices_; }

  /// set right_index
  void set_right_index(const std::vector<OrbitalIndex> &bra_idxs) {
    ket_indices_ = bra_idxs;
  }

  /// set left_index
  void set_left_index(const std::vector<OrbitalIndex> &ket_idxs) {
    bra_indices_ = ket_idxs;
  }

  /// set Operator
  void set_operator(const Operator &oper) { oper_ = oper; }

  /// set Notation
  void set_notation(const Notation &notation) {
    TA_USER_ASSERT(notation != Notation::Invalid, "invalid Notation")
    notation_ = notation;
  }

  /// set Operator type @sa Operator::Type
  void set_operator_type(const Operator::Type &oper_type) {
    oper_.set_type(oper_type);
  }

  /// Operator accessor
  const Operator &oper() const { return oper_; }

  /// Notation accessor
  const Notation &notation() const {
    TA_USER_ASSERT(notation_ != Notation::Invalid, "invalid Notation")
    return notation_;
  }

  /// @return true if this contains \c index
  bool has_index(const OrbitalIndex &index) const;

  /// @return true if it only contains AO indices
  bool is_ao() const;

  /// dimension of formula(2, 3 or 4)
  std::size_t rank() const;

  bool operator<(const Formula &other) const;
  bool operator==(const Formula &other) const;
  bool operator!=(const Formula &other) const { return !(*this == other); }

  /// converts this to a TA expression annotation
  /// @tparam Transformer a unary functor class
  /// @param transform_op used to transform index keys
  template <typename Transformer = detail::identity>
  std::string to_ta_expression(Transformer transform_op = detail::identity()) const;

 private:
  /// parse the index on one side
  std::vector<OrbitalIndex> check_orbital_index(std::wstring index_array);

 private:
  Operator oper_;
  Notation notation_;
  std::vector<OrbitalIndex> bra_indices_;
  std::vector<OrbitalIndex> ket_indices_;
};

template <typename Transformer>
std::string Formula::to_ta_expression(Transformer transform_op) const {
  std::string ta_expression;
  std::size_t rank = this->rank();
  std::size_t count = 0;

  // add left index
  for (const auto& index : bra_indices_) {
    std::string index_expression = transform_op(index.to_ta_expression());
    ta_expression.append(index_expression.begin(), index_expression.end());
    ++count;
    ta_expression.append(", ");
  }

  // add right index
  for (const auto& index : ket_indices_) {
    std::string index_expression = transform_op(index.to_ta_expression());
    ta_expression.append(index_expression.begin(), index_expression.end());
    ++count;
    if (count < rank) {
      ta_expression.append(", ");
    }
  }

  return ta_expression;
}

}

#endif  // MPQC_FORMULA_H
