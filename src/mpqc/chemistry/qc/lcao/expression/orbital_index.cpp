//
// Created by Chong Peng on 10/14/15.
//

#include <cstdlib>
#include <cwchar>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "mpqc/chemistry/qc/lcao/expression/greek_to_english_name.h"
#include "mpqc/chemistry/qc/lcao/expression/orbital_index.h"

namespace mpqc {
namespace lcao {

// set the range of key index to use
const wchar_t OrbitalIndex::frozen_occ_wchar[2] = {L'm', L'n'};
const wchar_t OrbitalIndex::occ_wchar[2] = {L'm', L'n'};
const wchar_t OrbitalIndex::corr_occ_wchar[2] = {L'i', L'l'};
const wchar_t OrbitalIndex::virt_wchar[2] = {L'a', L'd'};
const wchar_t OrbitalIndex::active_wchar[2] = {L'x', L'y'};
const wchar_t OrbitalIndex::any_wchar[2] = {L'p', L's'};
const wchar_t OrbitalIndex::othervirt_wchar[2] = {L'a', L'd'};
const wchar_t OrbitalIndex::allvirt_wchar[2] = {L'A', L'D'};
const wchar_t OrbitalIndex::allany_wchar[2] = {L'P', L'S'};
const wchar_t OrbitalIndex::obs_wchar[4] = {L'κ', L'λ', L'μ', L'ν'};
const wchar_t OrbitalIndex::vbs_wchar[4] = {L'Α', L'Β', L'Γ', L'Δ'};
const wchar_t OrbitalIndex::dfbs_wchar[4] = {L'Κ', L'Λ', L'Μ', L'Ν'};
const wchar_t OrbitalIndex::abs_wchar[4] = {L'α', L'β', L'γ', L'δ'};
const wchar_t OrbitalIndex::ribs_wchar[4] = {L'ρ', L'σ', L'τ', L'υ'};
const wchar_t OrbitalIndex::ubs_wchar[1] = {L'U'};

void OrbitalIndex::init(const std::wstring &string) {
  const wchar_t *letter = string.c_str();

  auto length = wcslen(letter);

  // single letter case, for example i, a , p

  const auto first = letter[0];
  if (length == 1) {
    index_ = wchar_to_index(first);
  }
  // letter with ' or letter with number, number is ignored
  // i1, a2, or A', a', P'
  else if (length == 2) {
    if (letter[1] == L'\'') {
      index_ = wchar_with_prime_to_index(first);
    } else if (isdigit(letter[1])) {
      index_ = wchar_to_index(first);
    } else {
      std::string error_message =
          "Wrong Key Index " + utility::to_string(string) + " ! \n";
      throw std::runtime_error(error_message);
    }
  }
  // letter with number and prime
  // A'1, a'1, P'12, i12
  else if (wcslen(letter) >= 3) {
    // with prime
    if (letter[1] == L'\'') {
      // check number
      for (auto i = 2ul; i < length; i++) {
        TA_ASSERT(isdigit(letter[i]));
      }

      index_ = wchar_with_prime_to_index(first);
    }
    // without prime
    else if (isdigit(letter[1])) {
      // check number
      for (auto i = 1ul; i < length; i++) {
        TA_ASSERT(isdigit(letter[i]));
      }

      index_ = wchar_to_index(first);
    } else {
      std::string error_message =
          "Wrong Key Index " + utility::to_string(string) + " ! \n";
      throw std::runtime_error(error_message);
    }
  }
}

bool OrbitalIndex::operator==(const OrbitalIndex &other) const {
  return (this->index_ == other.index_) && (this->spin_ == other.spin_);
}

bool OrbitalIndex::operator!=(const OrbitalIndex &other) const {
  return (this->index_ != other.index_) || (this->spin_ != other.spin_);
}

bool OrbitalIndex::operator<(const OrbitalIndex &other) const {
  return other.index() == this->index() ? this->spin() < other.spin()
                                        : this->index() < other.index();
}

bool OrbitalIndex::operator>(const OrbitalIndex &other) const {
  return !(*this < other) && !(*this == other);
}
bool OrbitalIndex::same(const OrbitalIndex &other) const {
  return (index_ == other.index()) && (this->spin_ == other.spin_) &&
         (name_ == other.name());
}

bool OrbitalIndex::is_ao() const {
  int index = static_cast<int>(index_);
  return index < 0;
}

bool OrbitalIndex::is_lcao() const {
  int index = static_cast<int>(index_);
  return index > 0;
}

OrbitalIndex::Type OrbitalIndex::wchar_to_index(const wchar_t first) {
  if (first >= occ_wchar[0] && first <= occ_wchar[1]) {
    return Type::occ;
  } else if (first >= corr_occ_wchar[0] && first <= corr_occ_wchar[1]) {
    return Type::corr_occ;
  } else if (first >= active_wchar[0] && first <= active_wchar[1]) {
    return Type::active;
  } else if (first >= virt_wchar[0] && first <= virt_wchar[1]) {
    return Type::virt;
  } else if (first >= any_wchar[0] && first <= any_wchar[1]) {
    return Type::any;
  } else if (first >= obs_wchar[0] && first <= obs_wchar[3]) {
    return Type::obs;
  } else if (first >= vbs_wchar[0] && first <= vbs_wchar[3]) {
    return Type::vbs;
  } else if (first >= abs_wchar[0] && first <= abs_wchar[3]) {
    return Type::abs;
  } else if (first >= dfbs_wchar[0] && first <= dfbs_wchar[3]) {
    return Type::dfbs;
  } else if (first >= ribs_wchar[0] && first <= ribs_wchar[3]) {
    return Type::ribs;
  } else if (first == ubs_wchar[0]) {
    return Type::ubs;
  } else {
    std::string error_message = "Wrong Key Index " +
                                utility::to_string(std::wstring(1, first)) +
                                " ! \n";
    throw std::runtime_error(error_message);
    return Type();
  }
}

OrbitalIndex::Type OrbitalIndex::wchar_with_prime_to_index(
    const wchar_t first) {
  if (first >= othervirt_wchar[0] && first <= othervirt_wchar[1]) {
    return Type::othervirt;
  } else if (first >= frozen_occ_wchar[0] && first <= frozen_occ_wchar[1]) {
    return Type::frozen_occ;
  } else if (first >= allvirt_wchar[0] && first <= allvirt_wchar[1]) {
    return Type::allvirt;
  } else if (first >= allany_wchar[0] && first <= allany_wchar[1]) {
    return Type::allany;
  } else {
    std::string error_message =
        "Wrong Key Index " + utility::to_string(std::wstring(1, first)) + " !";
    throw std::runtime_error(error_message);
    return Type();
  }
}

bool OrbitalIndex::is_mo_in_obs() const {
  int index = static_cast<int>(index_);
  return (index > 0) && (index <= 9);
}

bool OrbitalIndex::is_mo_in_abs() const { return index_ == Type::othervirt; }

bool OrbitalIndex::is_mo_in_ribs() const {
  int index = static_cast<int>(index_);
  return index >= 15;
}

OrbitalIndex OrbitalIndex::mo_to_ao() const {
  if (this->is_ao()) {
    return OrbitalIndex(*this);
  }

  std::wstring new_string;

  if (this->is_mo_in_abs()) {
    new_string.push_back(ribs_wchar[0]);
  } else if (this->is_mo_in_obs()) {
    new_string.push_back(obs_wchar[0]);
  } else if (this->is_mo_in_ribs()) {
    new_string.push_back(ribs_wchar[0]);
  }

  if (name_[1] == L'\'') {
    new_string.insert(new_string.end(), name_.begin() + 2, name_.end());
  } else {
    new_string.insert(new_string.end(), name_.begin() + 1, name_.end());
  }

  return OrbitalIndex(new_string);
}

std::string OrbitalIndex::to_ta_expression() const {
  std::string ta_expression;
  int length;
  const wchar_t *pt;
  pt = name_.c_str();
  char buffer[MB_CUR_MAX];

  while (*pt) {
    length = std::wctomb(buffer, *pt);
    TA_USER_ASSERT(
        length != -1,
        "OrbitalIndex::to_ta_expression encountered invalid character");
    // multiple byte, convert to english name
    if (length > 1) {
      auto pos = greek_to_english_name.find(*pt);
      if (pos == greek_to_english_name.end()) {
        throw std::runtime_error(
            "Couldn't Find The English Name of Greek Letter");
      }
      std::string convert_name = pos->second;
      ta_expression.append(convert_name.begin(), convert_name.end());
    }
    // single byte, do nothing
    else if (length == 1) {
      ta_expression.push_back(*buffer);
    }
    ++pt;
  }
  return ta_expression;
}

}  // namespace lcao
}  // namespace mpqc
