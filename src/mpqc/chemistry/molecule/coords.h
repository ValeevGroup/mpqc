#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_COORDS_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_COORDS_H_

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/util/misc/exenv.h"

namespace mpqc {

/** The MolecularCoordinates abstract class describes the coordinate system used
to describe a molecule.  It is used to convert a molecule's Cartesian
coordinates to and from this coordinate system. */
class MolecularCoordinates : virtual public DescribedClass {
 public:
  MolecularCoordinates(const std::shared_ptr<Molecule>& mol);
  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object to be queried
   *
   * The KeyVal object will be query the following keywords:
   * | Keyword | Type | Default| Description |
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

  /// clones this object
  /// @return a clone of this object
  virtual std::shared_ptr<MolecularCoordinates> clone() const = 0;

  /// prints this object to a std::ostream
  /// @params os the std::ostream object to which the output will be directed
  virtual void print(std::ostream& os = ExEnv::out0()) const = 0;

 protected:
  void update_molecule(const std::vector<Atom>& atoms) {
    molecule()->update(atoms);
  }

 private:
  std::shared_ptr<Molecule> molecule_;

  template <std::size_t N>
  friend void increment(MolecularCoordinates*, std::array<size_t, N>,
                        std::array<double, N>);

  virtual void displace(size_t ncoords, size_t* coords,
                        double* displacements) = 0;
};

template <std::size_t N>
void increment(MolecularCoordinates* coords, std::array<size_t, N> coord_idxs,
               std::array<double, N> step) {
  coords->displace(N, &coord_idxs[0], &step[0]);
}

template <std::size_t N>
void decrement(MolecularCoordinates* coords, std::array<size_t, N> coord_idxs,
               std::array<double, N> step) {
  std::array<double, N> inverse_step;
  for (size_t i = 0; i != N; ++i) inverse_step[i] = -step[i];
  increment(coords, coord_idxs, inverse_step);
}

std::ostream& operator<<(std::ostream& os, const MolecularCoordinates& coord);

/** The CartMolecularCoordinates class represents Cartesian coordinates of a
 * Molecule. */
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

  std::shared_ptr<MolecularCoordinates> clone() const override;

  void displace(size_t ncoords, size_t* coords, double* displacements) override;

  /// prints this object to a std::ostream
  /// @params os the std::ostream object to which the output will be directed
  void print(std::ostream& os = ExEnv::out0()) const override;
};

}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_MOLECULE_COORDS_H_ */
