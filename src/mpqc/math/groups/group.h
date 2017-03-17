/*
 * group.h
 *
 *  Created on: Feb 2, 2017
 *      Author: evaleev
 */

#ifndef SRC_MPQC_MATH_GROUPS_GROUP_H_
#define SRC_MPQC_MATH_GROUPS_GROUP_H_

#include <memory>
#include <vector>

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace math {

/// @brief Group is an abstract discrete group

/// Irreps can be indexed by arbitrary type (depends on the particular group,
/// e.g. a scalar
/// in ordinary point groups, but a vector of 3 reals in discretized space
/// groups). However,
/// the irreps are organized into a sequence, indexed by ordinal (scalar)
/// indices.
/// The range of ordinal indices is [0,nirreps) .
class Group : public std::enable_shared_from_this<Group> {
 public:
  using ordinal_type = std::size_t;

  /// an irreducible representation of Group
  class Irrep {
   public:
    virtual ~Irrep() { }

    /// @return true if this is a trivial representation
    virtual bool is_trivial() const = 0;

    /// @return shared_ptr to the Group object that defines this
    const std::weak_ptr<const Group>& group() const;

    /// tensor product of 2 Irrep objects = a "linear combination" of Irrep objects
    virtual std::vector<std::pair<unsigned int, std::shared_ptr<const Irrep>>>
    tensor_product(std::shared_ptr<const Irrep> other) const = 0;

   protected:
    Irrep(std::weak_ptr<const Group> grp);
    std::weak_ptr<const Group> group_;
  };

  /// a table of irreducible representations for Group
  class IrrepTable {
   public:
    virtual ~IrrepTable() { }

    /// @return shared_ptr to the Group object that defines this
    const std::weak_ptr<const Group>& group() const;

    /// @return the number of irreps
    ordinal_type size() const;

    /// @param irrep_ordinal the irrep ordinal
    /// @return the irrep object
    virtual std::shared_ptr<const Irrep> make_irrep(
        ordinal_type irrep_ordinal) const = 0;

    /// inverse of Group::irrep_index, converts irrep_idx to the irrep ordinal
    /// @param irrep_index irrep index, opaqued into a boost::any
    /// @return the irrep ordinal
    virtual ordinal_type irrep_ordinal(
        std::shared_ptr<const Irrep> irrep) const = 0;

   protected:
    IrrepTable(std::weak_ptr<const Group> grp);
    std::weak_ptr<const Group> group_;
  };

  virtual ~Group() { }

  /// @return the order of the group
  virtual ordinal_type order() const = 0;

  /// @return the IrrepTable object
  virtual std::shared_ptr<const IrrepTable> irrep_table() const = 0;
};

namespace groups {

/// @brief Z1 is the trivial group.
class Z1 : public Group {
 public:
  Z1() = default;
  ~Z1() = default;

  ordinal_type order() const override;

  std::shared_ptr<const IrrepTable> irrep_table() const override;

  class Irrep : public Group::Irrep {
   public:
    Irrep(std::shared_ptr<const Group> grp = std::make_shared<const Z1>());
    ~Irrep() = default;

    bool is_trivial() const override;

    std::vector<std::pair<unsigned int, std::shared_ptr<const Group::Irrep>>>
    tensor_product(std::shared_ptr<const Group::Irrep> other) const override;
  };

  class IrrepTable : public Group::IrrepTable {
   public:
    IrrepTable(std::shared_ptr<const Group> grp);
    ~IrrepTable() = default;

    std::shared_ptr<const Group::Irrep> make_irrep(
        ordinal_type irrep_ordinal) const override;
    ordinal_type irrep_ordinal(
        std::shared_ptr<const Group::Irrep> irrep) const override;
  };
};

/// @brief SupercellTranslationGroup is a discrete group of translations of a
/// lattice

/// SupercellTranslationGroup is a 3-d space group of lattice translations
/// modulo
/// a supercell composed of \f$ \mathbf{n} = \{n_a, n_b, n_c\} \f$ primitive
/// unit
/// cells of size \f$ \mathbf{L} = \{ L_a, L_b, L_c \} \f$ .
class SupercellTranslationGroup : public Group {
 public:
  /// @param L primitive unit cell dimensions, \f$ \mathbf{L} = \{ L_a, L_b, L_c
  /// \} \f$
  /// @param n supercell size, \f$ \mathbf{n} = \{n_a, n_b, n_c\} \f$
  SupercellTranslationGroup(Vector3d L, Vector3i n);

  ordinal_type order() const override;

  std::shared_ptr<const IrrepTable> irrep_table() const override;

  /// Irrep of SupercellTranslationGroup is labeled by wavenumbers
  /// \f$ \{k_a, k_b, k_c\} = \{ \frac{\pi (2 K_a - (n_a-1))}{L_a (n_a-1)},
  /// \frac{\pi (2 K_b - (n_b-1))}{L_b (n_b-1)}, \frac{\pi (2 K_c -
  /// (n_c-1))}{L_c (n_c-1)} \}
  /// \f$ where \f$ K_a = [0..n_a) \f$, etc. Internally, integers \f$ \{ K_a,
  /// K_b, K_c \} \f$
  /// are stored, and wave numbers are computed as needed.
  class Irrep : public Group::Irrep {
   public:
    /// constructs an Irrep, \c K is renormalized modulo \c n
    Irrep(std::shared_ptr<const Group> grp, Vector3i K);
    ~Irrep() = default;

    bool is_trivial() const override;

    std::vector<std::pair<unsigned int, std::shared_ptr<const Group::Irrep>>>
    tensor_product(std::shared_ptr<const Group::Irrep> other) const override;

    const Vector3i& K() const;

   private:
    Vector3i K_;

    /// converts weak_ptr to shared_ptr
    std::shared_ptr<const SupercellTranslationGroup> grp_() const;
  };

  class IrrepTable : public Group::IrrepTable {
   public:
    IrrepTable(std::shared_ptr<const Group>);
    ~IrrepTable() = default;
    /// converts an ordinal to the corresponding Irrep
    std::shared_ptr<const Group::Irrep> make_irrep(
        ordinal_type irrep_ordinal) const override;

    /// computes the ordinal correspondign to an Irrep
    ordinal_type irrep_ordinal(
        std::shared_ptr<const Group::Irrep> irrep) const override;

   private:
    /// converts weak_ptr to shared_ptr
    std::shared_ptr<const SupercellTranslationGroup> grp_() const;
  };

  /// @return (real space) unit cell dimensions
  const Vector3d& L() const;
  /// @return supercell size
  const Vector3i& n() const;

 private:
  Vector3d L_;  //!< (real space) unit cell dimensions
  Vector3i n_;  //!< supercell size
};

}  // namespace group
}  // namespace mpqc
}  // namespace mpqc

#endif /* SRC_MPQC_MATH_GROUPS_GROUP_H_ */
