#pragma once
#ifndef MPQC_MOLECULE_ATOM_H
#define MPQC_MOLECULE_ATOM_H

#include <string>
#include <iosfwd>

#include "./molecule_fwd.h"

namespace mpqc {
namespace molecule {

/*! \ingroup Molecule @{ */

/*! \brief A class which holds the basic information for an atom
 *
 *  Atom has a position in Bohr, an atomic number, and a mass.  By default the
 *  position is an Eigen::Vector3d.
 */
class Atom {
  private:
    Vec3D center_ = {0, 0, 0};
    int64_t atomic_number_ = 0;
    double mass_ = 0;

  public:
    Atom() = default;
    Atom(const Atom &atom) = default;
    Atom &operator=(const Atom &atom) = default;
    Atom(Atom &&atom) = default;
    Atom &operator=(Atom &&atom) = default;

    /*! \brief Constructs the atom with the center, mass and charge provided.
     */
    Atom(Vec3D const &center, double mass, int64_t Z)
            : center_(center), atomic_number_(Z), mass_(mass) {}

    /*! Returns the location of the atom in Bohr. */
    Vec3D const &center() const { return center_; }

    /*! Returns the charge of the atom in atomic units. */
    int64_t charge() const { return atomic_number_; }

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
inline Vec3D const &center(Atom const &a) { return a.center(); }

/// Returns the mass of the atom
inline double mass(Atom const &a) { return a.mass(); }

/// Returns the nuclear charge of the atom.
inline int64_t charge(Atom const &a) { return a.charge(); }

/*! \brief Returns the center of the atom.
 *
 * Center of mass is part of the clustering interface.
 */
inline Vec3D const &center_of_mass(Atom const &a) { return a.center(); }

/*! @} */

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_ATOM_H
