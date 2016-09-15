//
// Created by Chong Peng on 10/31/15.
//

#include "operator.h"

#include <boost/algorithm/string.hpp>
#include <memory>
#include <mpqc/util/misc/string.h>
#include <string>

namespace mpqc {

const std::map<Operator::Type, std::wstring> Operator::oper_to_string = {
    {Type::Overlap, L""},      {Type::Kinetic, L"T"},
    {Type::Nuclear, L"V"},     {Type::Core, L"H"},
    {Type::Coulomb, L"G"},     {Type::cGTG, L"R"},
    {Type::cGTG2, L"R2"},      {Type::cGTGCoulomb, L"GR"},
    {Type::DelcGTG2, L"dR2"},  {Type::J, L"J"},
    {Type::hJ, L"hJ"},         {Type::K, L"K"},
    {Type::KAlpha, L"K(α)"},   {Type::KBeta, L"K(β)"},
    {Type::Fock, L"F"},        {Type::FockAlpha, L"F(α)"},
    {Type::FockBeta, L"F(β)"}, {Type::Identity, L"I"}};

const std::map<Operator::Option, std::wstring> Operator::option_to_string = {
    {Operator::Option::DensityFitting, L"df"},
    {Operator::Option::Inverse, L"inv"},
    {Operator::Option::InverseSquareRoot, L"inv_sqr"}};

bool Operator::is_onebody() const {
  return Type::__first_1body_operator <= type_ &&
         type_ <= Type::__last_1body_operator;
}

bool Operator::is_twobody() const {
  return Type::__first_2body_operator <= type_ &&
         type_ <= Type::__last_2body_operator;
}

bool Operator::is_fock() const {
  return (type_ == Type::Fock || type_ == Type::FockAlpha ||
          type_ == Type::FockBeta);
}

bool Operator::is_jk() const {
  return (type_ == Type::J || type_ == Type::K || type_ == Type::KAlpha ||
          type_ == Type::KBeta);
}

bool Operator::is_r12() const {
  return (type_ == Type::cGTG2 || type_ == Type::cGTG ||
          type_ == Type::cGTGCoulomb || type_ == Type::DelcGTG2);
}

bool Operator::has_option(Operator::Option op) const {
  auto df = std::find(option_.cbegin(), option_.cend(), op);
  return (df != option_.cend());
}

Operator::Operator(std::wstring oper, std::wstring opt) {
  // parse operation
  if (oper.empty()) {
    throw std::runtime_error(utility::to_string(oper) + " Empty Operation! \n");
  }

  boost::trim(oper);

  auto iter = std::find_if(
      begin(oper_to_string), end(oper_to_string),
      [=](const std::pair<Operator::Type, std::wstring> item) -> bool {
        return item.second == oper;
      });

  if (iter == oper_to_string.end())
    throw std::runtime_error(utility::to_string(oper) +
                             " Invalid Operation! \n");
  else
    type_ = iter->first;

  // parse option

  // split string by , and space
  if (!opt.empty()) {
    std::vector<std::wstring> split_option;
    boost::split(split_option, opt, boost::is_any_of(L", "),
                 boost::token_compress_on);

    std::vector<Operator::Option> result;
    for (const auto& op : split_option) {
      auto iter = std::find_if(
          begin(option_to_string), end(option_to_string),
          [=](const std::pair<Operator::Option, std::wstring> item) -> bool {
            return item.second == op;
          });
      if (iter == option_to_string.end()) {
        throw std::runtime_error(utility::to_string(op) +
                                 " Invalid Option! \n");
      } else {
        result.push_back(iter->first);
      }
    }
    std::sort(result.begin(), result.end());
    option_ = result;
  }
}

bool Operator::operator==(const Operator& other) const {
  bool same_operation = (this->type_ == other.type_);

  bool same_option = false;

  if (option_.size() == other.option_.size()) {
    same_option =
        std::equal(option_.begin(), option_.end(), other.option_.begin());
  }

  return same_operation && same_option;
}

bool Operator::operator<(const Operator& other) const {
  return other.type() == this->type() ? this->option() < other.option()
                                      : this->type() < other.type();
}

const std::wstring Operator::oper_string() const {
  const auto result = oper_to_string.find(type_);
  return result->second;
}

const std::wstring Operator::option_string() const {
  std::wstring result;
  if (option_.empty()) {
    return result;
  }
  for (const auto& option : option_) {
    result += option_to_string.find(option)->second + L",";
  }
  result = L"[" + result;
  result.back() = L']';
  return result;
}
}  // namespace mpqc
