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

  /// check if formula has index in left_index and right_index
  bool has_index(const OrbitalIndex &index) const;

  /// dimension of formula(2, 3 or 4)
  std::size_t rank() const;

  bool operator<(const Formula &other) const;
  bool operator==(const Formula &other) const;
  bool operator!=(const Formula &other) const { return !(*this == other); }

  /// convert to TA expression string format
  std::string to_ta_expression() const;

 private:
  /// parse the index on one side
  std::vector<OrbitalIndex> check_orbital_index(std::wstring index_array);

 private:
  Operator oper_;
  Notation notation_;
  std::vector<OrbitalIndex> bra_indices_;
  std::vector<OrbitalIndex> ket_indices_;
};

}

#endif  // MPQC_FORMULA_H
