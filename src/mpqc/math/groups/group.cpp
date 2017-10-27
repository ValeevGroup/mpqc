
#include "mpqc/math/groups/group.h"

#include <cassert>
#include <memory>

#include "mpqc/util/core/exception.h"

namespace mpqc {
namespace math {

Group::Irrep::Irrep(std::weak_ptr<const Group> grp) : group_(grp) {}

const std::weak_ptr<const Group>& Group::Irrep::group() const { return group_; }

Group::ordinal_type Group::IrrepTable::size() const {
  if (auto grp = group_.lock()) {
    return grp->order();
  } else
    throw ProgrammingError("Group::IrrepTable: bad group ptr", __FILE__,
                           __LINE__);
}

Group::IrrepTable::IrrepTable(std::weak_ptr<const Group> grp) : group_(grp) {}

const std::weak_ptr<const Group>& Group::IrrepTable::group() const {
  return group_;
}

namespace groups {

//////////////////////////////////////////////////////

Group::ordinal_type Z1::order() const { return 1; }

std::shared_ptr<const Group::IrrepTable> Z1::irrep_table() const {
  return std::make_shared<Z1::IrrepTable>(this->shared_from_this());
}

Z1::Irrep::Irrep(std::shared_ptr<const Group> grp) : Group::Irrep(grp) {}

bool Z1::Irrep::is_trivial() const { return true; }

std::vector<std::pair<unsigned int, std::shared_ptr<const Group::Irrep>>>
Z1::Irrep::tensor_product(std::shared_ptr<const Group::Irrep> other) const {
  return {std::make_pair(1, std::make_shared<Z1::Irrep>())};
}

Z1::IrrepTable::IrrepTable(std::shared_ptr<const Group> grp) : Group::IrrepTable(grp) {
}

std::shared_ptr<const Group::Irrep> Z1::IrrepTable::make_irrep(
    ordinal_type irrep_ordinal) const {
  assert(irrep_ordinal == 0);
  return std::make_shared<const Z1::Irrep>();
}

Group::ordinal_type Z1::IrrepTable::irrep_ordinal(
    std::shared_ptr<const Group::Irrep> irrep) const {
  assert(std::dynamic_pointer_cast<const Z1::Irrep>(irrep) != nullptr);
  return 0;
}

//////////////////////////////////////////////////////

SupercellTranslationGroup::SupercellTranslationGroup(Vector3d L, Vector3i n)
    : L_(L), n_(n) {}

Group::ordinal_type SupercellTranslationGroup::order() const {
  return n_(0) * n_(1) * n_(2);
}

const Vector3d& SupercellTranslationGroup::L() const { return L_; }

const Vector3i& SupercellTranslationGroup::n() const { return n_; }

std::shared_ptr<const Group::IrrepTable>
SupercellTranslationGroup::irrep_table() const {
  return std::make_shared<const SupercellTranslationGroup::IrrepTable>(
      this->shared_from_this());
}

SupercellTranslationGroup::Irrep::Irrep(std::shared_ptr<const Group> grp,
                                        Vector3i K)
    : Group::Irrep(grp), K_(K) {
  auto n = grp_()->n();
  for (int i = 0; i != 3; ++i) K_(i) = K_(i) % n(i);
}

bool SupercellTranslationGroup::Irrep::is_trivial() const {
  auto n = grp_()->n();
  // k == 0 ~ 2K - (n-1) == 0
  return (2 * K_ - n == Vector3i{-1, -1, -1});
}

std::vector<std::pair<unsigned int, std::shared_ptr<const Group::Irrep>>>
SupercellTranslationGroup::Irrep::tensor_product(
    std::shared_ptr<const Group::Irrep> other) const {
  auto other_ptr =
      std::static_pointer_cast<const SupercellTranslationGroup::Irrep>(other);
  assert(this->grp_() == other_ptr->grp_());
  return {std::make_pair(
      1, std::make_shared<const Irrep>(grp_(), K() + other_ptr->K()))};
}

const Vector3i& SupercellTranslationGroup::Irrep::K() const { return K_; }

std::shared_ptr<const SupercellTranslationGroup>
SupercellTranslationGroup::Irrep::grp_() const {
  if (auto grp_shr = this->group().lock())
    return std::static_pointer_cast<const SupercellTranslationGroup>(grp_shr);
  else
    throw ProgrammingError("SupercellTranslationGroup::Irrep: bad group ptr",
                           __FILE__, __LINE__);
}

SupercellTranslationGroup::IrrepTable::IrrepTable(
    std::shared_ptr<const Group> grp)
    : Group::IrrepTable(grp) {}

std::shared_ptr<const SupercellTranslationGroup>
SupercellTranslationGroup::IrrepTable::grp_() const {
  if (auto grp_shr = this->group().lock())
    return std::static_pointer_cast<const SupercellTranslationGroup>(grp_shr);
  else
    throw ProgrammingError("SupercellTranslationGroup::Irrep: bad group ptr",
                           __FILE__, __LINE__);
}

Group::ordinal_type SupercellTranslationGroup::IrrepTable::irrep_ordinal(
    std::shared_ptr<const Group::Irrep> irrep) const {
  auto irrep_ptr =
      std::static_pointer_cast<const SupercellTranslationGroup::Irrep>(irrep);
  auto K = irrep_ptr->K();
  auto n = grp_()->n();
  return (K(0) * n(1) + K(1)) * n(2) + K(2);
}

std::shared_ptr<const Group::Irrep>
SupercellTranslationGroup::IrrepTable::make_irrep(
    ordinal_type irrep_ordinal) const {
  assert(irrep_ordinal < size());
  auto n = grp_()->n();
  auto n12 = n(1) * n(2);
  Vector3i K;
  K(0) = irrep_ordinal / n12;
  auto K1_n2_plus_K2 = irrep_ordinal % n12;
  K(1) = K1_n2_plus_K2 / n(2);
  K(2) = K1_n2_plus_K2 % n(2);
  return std::make_shared<SupercellTranslationGroup::Irrep>(grp_(), K);
}

}  // namespace groups
}  // namespace math
}  // namespace mpqc
