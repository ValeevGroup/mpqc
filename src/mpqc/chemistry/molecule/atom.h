
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_H_

#include <cmath>
#include <iosfwd>
#include <string>

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {

/*! \ingroup Molecule @{ */

/*! \brief A class which holds the basic information for an atom
 *
 *  Atom has a position in Bohr, an atomic number, and a mass.  By default the
 *  position is an Eigen::Vector3d.
 */
class Atom {
 private:
  Vector3d center_ = {0, 0, 0};
  int64_t atomic_number_ = 0;
  double charge_ = 0.0;
  double mass_ = 0;

 public:
  Atom() = default;
  Atom(const Atom &atom) = default;
  Atom &operator=(const Atom &atom) = default;
  Atom(Atom &&atom) = default;
  Atom &operator=(Atom &&atom) = default;

  // clang-format off
  /// \brief The KeyVal constructor
  /// \param kv The KeyVal object. The following keywords will be queried:
  ///  | Keyword | Type | Default| Description |
  ///  |---------|------|--------|-------------|
  ///  |\c element|string|none|the element symbol (e.g., H, He, etc.).|
  ///  |\c xyz|array|none|a 3-element array of xyz coordinates.|
  ///  |\c mass|real|mass of the most common isotope|the mass in amu.|
  ///  |\c charge|real|the atomic number|nuclear charge|
  // clang-format on
  Atom(const KeyVal &kv);

  /*! \brief Constructs the atom with the center, mass and atomix number
   * provided.
   */
  Atom(Vector3d const &center, double mass, int64_t Z, double charge = NAN)
      : center_(center),
        atomic_number_(Z),
        charge_(std::isnan(charge) ? static_cast<double>(Z) : charge),
        mass_(mass)
  {}

  /*! @return the Cartesian coordinates (in atomic units). */
  Vector3d const &center() const { return center_; }

  /*! @return the atomic number. */
  int64_t atomic_number() const { return atomic_number_; }

  /*! @return the nuclearcharge */
  double charge() const { return charge_; }

  /*! Returns the mass of the atom in atomic units. */
  double mass() const { return mass_; }

  /*! \brief Returns the atom in xyz format.
   *
   * By default the function assumes the units are in Bohr, so it will
   * convert to ang, but by passing false to the function it will
   * not convert the position.
   */
  std::string xyz_string(bool convert_to_angstroms = true) const;
};

//
// External interface functions
//

std::ostream &operator<<(std::ostream &, Atom const &);

/// Returns the center of the atom
inline Vector3d const &center(Atom const &a) { return a.center(); }

/// Returns the mass of the atom
inline double mass(Atom const &a) { return a.mass(); }

/// @returns the atomic number of the Atom object
inline int64_t total_atomic_number(Atom const &a) { return a.atomic_number(); }

/// # of atoms in an atom, always returns 1
inline size_t natoms(Atom const &a) { return 1; }

/*! \brief Returns the center of the atom.
 *
 * Center of mass is part of the clustering interface.
 */
inline Vector3d const &center_of_mass(Atom const &a) { return a.center(); }

/*! @} */

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_H_
