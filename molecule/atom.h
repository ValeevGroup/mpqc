#pragma once
#ifndef TCC_MOLECULE_ATOM_H
#define TCC_MOLECULE_ATOM_H

#include <array>
#include "molecule_fwd.h"

namespace tcc {
namespace molecule {

/*! \class Atom
 *  \brief A class which holds the basic information for an atom
 *
 *  Atom has a position, an atomic number, and a mass.  By default the 
 *  position is an Eigen::Vector3d.
 */
class Atom {
  public:
    Atom() = default;
    Atom(const Atom &atom) = default;
    Atom &operator=(const Atom &atom) = default;
    Atom(Atom &&atom) = default;
    Atom &operator=(Atom &&atom) = default;

    Atom(position_t center, double mass, int Z)
        : center_(std::move(center)), atomic_number_(Z), mass_(mass) {}

    /*! Returns the location of the atom in Bohr. */
    position_t center() const { return center_; }
    /*! Returns the charge of the atom in atomic units. */
    int charge() const { return atomic_number_; }
    /*! Returns the mass of the atom in atomic units. */
    double mass() const { return mass_; }

  private:
    position_t center_ = {0, 0, 0};
    int atomic_number_ = 0;
    double mass_ = 0;
};

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_ATOM_H
