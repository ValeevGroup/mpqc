//
// Created by Chong Peng on 10/31/15.
//

#include "operator.h"

#include "greek_to_english_name.h"
#include "mpqc/util/misc/string.h"
#include <boost/algorithm/string.hpp>

#include <codecvt>
#include <locale>
#include <memory>
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
    {Type::FockBeta, L"F(β)"}, {Type::Cadf, L"Cadf"},
    {Type::Identity, L"I"}};

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

bool Operator::is_jk() const { return (type_ == Type::J || type_ == Type::K || type_ == Type::KAlpha ||
          type_ == Type::KBeta);
}

bool Operator::is_r12() const {
  return (type_ == Type::cGTG2 || type_ == Type::cGTG ||
          type_ == Type::cGTGCoulomb || type_ == Type::DelcGTG2);
}

Operator::Operator(std::wstring oper) {
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
}

bool Operator::operator==(const Operator& other) const {
  return (this->type_ == other.type_);
}

bool Operator::operator<(const Operator& other) const {
  return (this->type() < other.type());
}

const std::wstring Operator::oper_string() const {
  const auto result = oper_to_string.find(type_);
  return result->second;
}

}  // namespace mpqc
