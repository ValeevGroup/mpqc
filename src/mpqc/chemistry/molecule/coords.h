#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_COORDS_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_COORDS_H_

#include "mpqc/chemistry/molecule/molecule.h"

namespace mpqc {

/** The MolecularCoordinates abstract class describes the coordinate system used
to describe a molecule.  It is used to convert a molecule's Cartesian
coordinates to and from this coordinate system. */
class MolecularCoordinates : public DescribedClass {
 public:
  MolecularCoordinates(const std::shared_ptr<Molecule>& mol);
  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object to be queried
   *
   * The KeyVal object will be query the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | molecule | Molecule | none | the molecule  |
   */
  // clang-format on
  MolecularCoordinates(const KeyVal& kv);

  virtual ~MolecularCoordinates();

  /// Reports the number of coordinates
  /// @return the number of coordinates
  virtual size_t size() const = 0;

  /// Reports the number of constrained coordinates. The default is zero.
  /// @return the number of constrained coordinates
  virtual size_t nconstrained() const;

  /// @return pointer to the molecule
  const std::shared_ptr<Molecule> molecule() const { return molecule_; }

 private:
  std::shared_ptr<Molecule> molecule_;
};

/** The CartMolecularCoordinates class represents Cartesian coordinates of a Molecule. */
class CartMolecularCoordinates : public MolecularCoordinates {
 public:
  CartMolecularCoordinates(const std::shared_ptr<Molecule>& mol);

  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object to be queried
   *
   * The KeyVal object will be query the keywords of class MolecularCoordinates.
   */
  // clang-format on
  CartMolecularCoordinates(const KeyVal& kv);

  virtual ~CartMolecularCoordinates();

  /// Reports the number of coordinates
  /// @return the number of coordinates
  size_t size() const override;
};

}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_MOLECULE_COORDS_H_ */
