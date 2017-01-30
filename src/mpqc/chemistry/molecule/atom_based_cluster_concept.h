
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_BASED_CLUSTER_CONCEPT_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_BASED_CLUSTER_CONCEPT_H_

#include "mpqc/chemistry/molecule/cluster_collapse.h"
#include "mpqc/chemistry/molecule/cluster_concept.h"
#include "mpqc/chemistry/molecule/molecule_fwd.h"

#include <iostream>
#include <memory>
#include <vector>

namespace mpqc {

/*!
 * \addtogroup ChemistryMolecule
 *
 * @{
 */

/*
 * ClusterConcept is the base class which defines the operations which
 * different clusterable types must have.
 */
class AtomBasedClusterConcept : public ClusterConcept {
 public:
  virtual ~AtomBasedClusterConcept() noexcept = default;
  virtual int64_t charge_() const = 0;
  virtual double mass_() const = 0;
  virtual Vector3d const &com_() const = 0;
  virtual std::vector<Atom> atoms_() const = 0;
  virtual size_t natoms_() const = 0;
  virtual void update_(const std::vector<Atom>& atoms, size_t& pos) = 0;
};

/*
 * AtomBasedClusterModel is a class which is basically used for type erasure
 */
template <typename T>
class AtomBasedClusterModel : public AtomBasedClusterConcept {
 private:
  T element_;

 public:
  AtomBasedClusterModel(T t) : element_(std::move(t)) {}
  AtomBasedClusterModel(const AtomBasedClusterModel &c) = default;
  AtomBasedClusterModel &operator=(AtomBasedClusterModel c) {
    element_ = std::move(c.element_);
    return *this;
  }

  AtomBasedClusterModel(AtomBasedClusterModel &&c) = default;
  AtomBasedClusterModel &operator=(AtomBasedClusterModel &&c) = default;

  AtomBasedClusterConcept *clone_() const override final {
    return new AtomBasedClusterModel(*this);
  }

  Vector3d const &center_() const override final { return center(element_); }
  Vector3d const &com_() const override final {
    return center_of_mass(element_);
  }

  int64_t charge_() const override final { return charge(element_); }
  double mass_() const override final { return mass(element_); }

  std::vector<Atom> atoms_() const override final {
    return collapse_to_atoms(element_);
  }

  void update_(const std::vector<Atom>& atoms, size_t& pos) override final {
    update(element_, atoms, pos);
  }

  size_t natoms_() const override final { return natoms(element_); }

  std::ostream &print_(std::ostream &os) const override final {
    os << element_;
    return os;
  }
};

/*!
 * \brief The AtomBasedClusterable is a class that holds any clusterable type
 * that is built up from atoms.
 *
 * AtomBasedClusterables must be able to return the total mass of
 * the clusterable as well as the total charge. Finally they need to be
 * collapsable to a vector of atoms.
 */
class AtomBasedClusterable {
 private:
  std::shared_ptr<AtomBasedClusterConcept> element_impl_;

 public:
  template <typename C>
  explicit AtomBasedClusterable(C c)
      : element_impl_(
            std::make_shared<AtomBasedClusterModel<C>>(std::move(c))) {}
  AtomBasedClusterable(AtomBasedClusterable const &c) = default;
  AtomBasedClusterable &operator=(AtomBasedClusterable const &c) = default;
  AtomBasedClusterable(AtomBasedClusterable &&c) = default;
  AtomBasedClusterable &operator=(AtomBasedClusterable &&c) = default;

  // Don't provide a way to access the center
  // Vector3d const &center() const { return element_impl_->center_(); }
  Vector3d const &com() const { return element_impl_->com_(); }

  /// Vector of atoms that make up the clusterable
  std::vector<Atom> atoms() const { return element_impl_->atoms_(); }

  double mass() const { return element_impl_->mass_(); }
  double charge() const { return element_impl_->charge_(); }
  double natoms() const { return element_impl_->natoms_(); }

  void update(const std::vector<Atom> &atoms, size_t &pos) const {
    return element_impl_->update_(atoms, pos);
  }

  std::ostream &print(std::ostream &os) const {
    return element_impl_->print_(os);
  }
};

inline double mass(AtomBasedClusterable const &ac) { return ac.mass(); }

inline double charge(AtomBasedClusterable const &ac) { return ac.charge(); }

inline size_t natoms(AtomBasedClusterable const &ac) { return ac.natoms(); }

inline Vector3d const &center(AtomBasedClusterable const &ac) {
  return ac.com();
}

inline Vector3d const &center_of_mass(AtomBasedClusterable const &ac) {
  return ac.com();
}

inline std::vector<Atom> collapse_to_atoms(AtomBasedClusterable const &ac) {
  return ac.atoms();
}

inline void update(AtomBasedClusterable &ac, const std::vector<Atom>& atoms, size_t& pos) {
  return ac.update(atoms, pos);
}

/*! @} */

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_BASED_CLUSTER_CONCEPT_H_
