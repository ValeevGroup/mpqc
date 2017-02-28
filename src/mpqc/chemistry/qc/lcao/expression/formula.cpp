//
// Created by Chong Peng on 10/15/15.
//

#include <algorithm>
#include <iostream>

#include "mpqc/chemistry/qc/lcao/expression/formula.h"

#include <TiledArray/error.h>
#include <boost/algorithm/string.hpp>

#include "mpqc/util/misc/exception.h"

namespace mpqc {

const std::map<Formula::Option, std::wstring> Formula::option_to_string = {
    {Formula::Option::DensityFitting, L"df"},
    {Formula::Option::Inverse, L"inv"},
    {Formula::Option::InverseSquareRoot, L"inv_sqr"}};

Formula::Formula(std::wstring string) {
  if (string.empty()) {
    throw std::runtime_error(utility::to_string(string) + " Empty Formula! \n");
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
        utility::to_string(string) +
        " Formula should start with < or ( and end with > or ). \n");
  }

  TA_ASSERT((*bra_symbol == L'<' && *ket_symbol == L'>') ||
            (*bra_symbol == L'(') && *ket_symbol == L')');

  if (*bra_symbol == L'<') {
    notation_ = Notation::Physical;
  } else {
    notation_ = Notation::Chemical;
  }

  std::wstring formula(bra_symbol + 1, ket_symbol);

  // find option between [ and ]
  auto option_left = std::find(string.cbegin(), string.cend(), L'[');
  auto option_right = std::find(string.cbegin(), string.cend(), L']');

  std::wstring opt;

  if ((option_left != string.cend()) && (option_right != string.cend())) {
    TA_ASSERT(option_left < option_right);
    opt = std::wstring(option_left + 1, option_right);
  }

  // parse options (separators = {",", " "})
  if (!opt.empty()) {
    std::vector<std::wstring> split_option;
    boost::split(split_option, opt, boost::is_any_of(L", "),
                 boost::token_compress_on);

    std::vector<Formula::Option> result;
    for (const auto& op : split_option) {
      // check if this is one of Options
      auto opt_iter = std::find_if(
          begin(option_to_string), end(option_to_string),
          [=](const std::pair<Formula::Option, std::wstring> item) -> bool {
            return item.second == op;
          });
      const auto has_opt = (opt_iter != option_to_string.end());

      // or math::PetiteList::Symmetry
      auto op_str = utility::to_string(op);
      auto symm_iter = std::find_if(
          begin(math::PetiteList::symmetry_to_string), end(math::PetiteList::symmetry_to_string),
          [=](const std::pair<math::PetiteList::Symmetry, std::string> item) -> bool {
            return item.second == op_str;
          });
      const auto has_symm = (symm_iter != math::PetiteList::symmetry_to_string.end());

      if (!has_opt && !has_symm)
        throw ProgrammingError((utility::to_string(op) +
                                 " Invalid Option! \n").c_str(), __FILE__, __LINE__);
      if (has_opt)
        result.push_back(opt_iter->first);
      if (has_symm)
        symm_ = symm_iter->first;
    }
    std::sort(result.begin(), result.end());
    options_ = result;
  }

  // split the formula by |
  std::vector<std::wstring> split_formula;
  boost::split(split_formula, formula, boost::is_any_of(L"|"),
               boost::token_compress_on);

  if (split_formula.size() == 2) {
    std::wstring left_formula = std::move(split_formula[0]);
    std::wstring right_formula = std::move(split_formula[1]);

    oper_ = Operator(L" ");
    bra_indices_ = check_orbital_index(left_formula);
    ket_indices_ = check_orbital_index(right_formula);

  } else if (split_formula.size() == 3) {
    std::wstring left_formula = std::move(split_formula[0]);
    std::wstring operation = std::move(split_formula[1]);
    std::wstring right_formula = std::move(split_formula[2]);

    oper_ = Operator(operation);
    bra_indices_ = check_orbital_index(left_formula);
    ket_indices_ = check_orbital_index(right_formula);

  } else {
    throw std::runtime_error(utility::to_string(string) +
                             " Formula in wrong length. \n");
  }

  // error detecting

  // one body operation
  if (oper_.is_onebody() && (bra_indices_.size() != 1) &&
      (ket_indices_.size() != 1)) {
    throw std::runtime_error(utility::to_string(string) +
                             "One body Operator with Wrong Index Size! \n");
  }

  // more than three index
  if ((bra_indices_.size() >= 3) || (ket_indices_.size() >= 3)) {
    throw std::runtime_error(utility::to_string(string) +
                             " Wrong Number of Index! \n");
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
  } else if (options_ != other.options_) {
    return options_ < other.options_;
  } else if (symm_ != other.symm_) {
    return symm_ < other.symm_;
  } else if ((this->rank() != 2) && (notation_ != other.notation())) {
    return notation_ < other.notation();
  } else if (bra_indices() != other.bra_indices()) {
    return bra_indices() < other.bra_indices();
  } else {
    return ket_indices() < other.ket_indices();
  }
}

bool Formula::operator==(const Formula& other) const {
  // special case
  if (this->rank() == 2) {
    return (oper_ == other.oper_) && (bra_indices_ == other.bra_indices_) &&
           (ket_indices_ == other.ket_indices_);
  } else {
    return (oper_ == other.oper_) && (options_ == other.options_) &&
           (symm_ == other.symm_) && (bra_indices_ == other.bra_indices_) &&
           (ket_indices_ == other.ket_indices_) &&
           (notation_ == other.notation_);
  }
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

  // append options
  auto option_string = [=]() -> std::wstring {
    std::wstring result;
    const auto has_options = !options_.empty();
    const auto has_symm = symm_ != math::PetiteList::Symmetry::e;
    if (!has_options && !has_symm) {
      return result;
    }
    result = L"[";
    for (const auto& option : options_) {
      result += option_to_string.find(option)->second + L",";
    }
    // include symmetry, if not trivial
    if (has_symm)
      result += utility::to_wstring(to_string(symm_)) + L",";
    result.back() = L']';
    return result;
  };

  result = result + option_string();

  return result;
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

bool Formula::is_ao() const {
  for (const auto& index : bra_indices_) {
    if (not index.is_ao()) {
      return false;
    }
  }
  for (const auto& index : ket_indices_) {
    if (not index.is_ao()) {
      return false;
    }
  }
  return true;
}

void Formula::set_option(Option op) {
  if (!has_option(op))
    options_.push_back(op);
}

bool Formula::has_option(Formula::Option op) const {
  auto df = std::find(options_.cbegin(), options_.cend(), op);
  return (df != options_.cend());
}

}  // namespace mpqc
