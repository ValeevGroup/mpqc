//
// Created by Chong Peng on 10/15/15.
//

#include <algorithm>
#include <iostream>

#include <mpqc/chemistry/qc/expression/formula.h>

#include <TiledArray/error.h>
#include <boost/algorithm/string.hpp>

namespace mpqc {

Formula::Formula(std::wstring string) {
  if (string.empty()) {
    throw std::runtime_error("Empty Formula!");
  }

  // detect the brackets <> or ()
  auto left_mark = [](const wchar_t letter) {
    return letter == L'<' || letter == L'(';
  };

  // find formula between < > or ( )
  auto bra_symbol = std::find_if(string.cbegin(), string.cend(), left_mark);

  auto ket_symbol = string.cbegin();

  if (*bra_symbol == L'<') {
    auto pos = string.find_last_of(L'>');
    ket_symbol += pos;
  } else {
    auto pos = string.find_last_of(L')');
    ket_symbol += pos;
  }

  TA_ASSERT(bra_symbol < ket_symbol);

  if (bra_symbol == string.cend() || ket_symbol == string.cend()) {
    throw std::runtime_error(
        "Formula should start with < or ( and end with > or ). \n");
  }

  TA_ASSERT((*bra_symbol == L'<' && *ket_symbol == L'>') ||
            (*bra_symbol == L'(') && *ket_symbol == L')');

  if (*bra_symbol == L'<') {
    notation_ = Notation::Physical;
  } else {
    notation_ = Notation::Chemical;
  }

  std::wstring formula(bra_symbol + 1, ket_symbol);

  //        std::wcout << formula_string;

  // find option between [ and ]
  auto option_left = std::find(string.cbegin(), string.cend(), L'[');
  auto option_right = std::find(string.cbegin(), string.cend(), L']');

  std::wstring option;

  if ((option_left != string.cend()) && (option_right != string.cend())) {
    TA_ASSERT(option_left < option_right);
    option = std::wstring(option_left + 1, option_right);
  }

  // split the formula by |
  std::vector<std::wstring> split_formula;
  boost::split(split_formula, formula, boost::is_any_of(L"|"),
               boost::token_compress_on);

  if (split_formula.size() == 2) {
    std::wstring left_formula = std::move(split_formula[0]);
    std::wstring right_formula = std::move(split_formula[1]);

    oper_ = Operator(L" ", option);
    bra_indices_ = check_orbital_index(left_formula);
    ket_indices_ = check_orbital_index(right_formula);

  } else if (split_formula.size() == 3) {
    std::wstring left_formula = std::move(split_formula[0]);
    std::wstring operation = std::move(split_formula[1]);
    std::wstring right_formula = std::move(split_formula[2]);

    oper_ = Operator(operation, option);
    bra_indices_ = check_orbital_index(left_formula);
    ket_indices_ = check_orbital_index(right_formula);

  } else {
    throw std::runtime_error("Formula in wrong length. \n");
  }

  // error detecting

  // one body operation
  if (oper_.is_onebody() && (bra_indices_.size() != 1) &&
      (ket_indices_.size() != 1)) {
    throw std::runtime_error("One body Operator with Wrong Index Size!");
  }

  // more than three index
  if ((bra_indices_.size() >= 3) || (ket_indices_.size() >= 3)) {
    throw std::runtime_error("Wrong Number of Index!");
  }
}

std::size_t Formula::rank() const {
  return (bra_indices_.size() + ket_indices_.size());
}

std::vector<OrbitalIndex> Formula::check_orbital_index(
    std::wstring index_array) {
  std::vector<std::wstring> split_index;
  boost::split(split_index, index_array, boost::is_any_of(L" "),
               boost::token_compress_on);

  if (split_index.empty()) {
    return std::vector<OrbitalIndex>();
  } else {
    std::vector<OrbitalIndex> result;

    for (auto index : split_index) {
      if (!index.empty()) {
        result.emplace_back(index);
      }
    }
    return result;
  }
}

bool Formula::operator<(const Formula& other) const {
  if (oper() != other.oper()) {
    return oper() < other.oper();
  } else if (notation_ != other.notation()) {
    return notation_ < other.notation();
  } else if (bra_indices() != other.bra_indices()) {
    return bra_indices() < other.bra_indices();
  } else {
    return ket_indices() < other.ket_indices();
  }
}

bool Formula::operator==(const Formula& other) const {
  return (oper_ == other.oper_) && (bra_indices_ == other.bra_indices_) &&
         (ket_indices_ == other.ket_indices_) && (notation_ == other.notation_);
}

std::wstring Formula::string() const {
  std::wstring left, right, result;
  for (const auto& index : bra_indices_) {
    left += index.name() + L" ";
  }

  for (const auto& index : ket_indices_) {
    right += L" " + index.name();
  }

  // add operation
  auto oper_str = oper_.oper_string();
  if (!oper_str.empty()) {
    result = left + L"|" + oper_.oper_string() + L"|" + right;
  } else {
    result = left + L"|" + right;
  }

  if (notation_ == Notation::Chemical) {
    result = L"( " + result + L" )";
  } else {
    result = L"< " + result + L" >";
  }
  // add option
  result = result + oper_.option_string();

  return result;
}

std::string Formula::to_ta_expression() const {
  std::string ta_expression;
  std::size_t rank = this->rank();
  std::size_t count = 0;

  // add left index
  for (const auto& index : bra_indices_) {
    std::string index_expression = index.to_ta_expression();
    ta_expression.append(index_expression.begin(), index_expression.end());
    ++count;
    ta_expression.append(", ");
  }

  // add right index
  for (const auto& index : ket_indices_) {
    std::string index_expression = index.to_ta_expression();
    ta_expression.append(index_expression.begin(), index_expression.end());
    ++count;
    if (count < rank) {
      ta_expression.append(", ");
    }
  }

  return ta_expression;
}

bool Formula::has_index(const OrbitalIndex& ind) const {
  for (auto& index : bra_indices_) {
    if (ind == index) {
      return true;
    }
  }

  for (auto& index : ket_indices_) {
    if (ind == index) {
      return true;
    }
  }

  return false;
}
}
