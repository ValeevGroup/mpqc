#pragma once
#ifndef MPQC_MOLECULE_MOLECULE_H
#define MPQC_MOLECULE_MOLECULE_H

#include "molecule_fwd.h"
#include <vector>

namespace mpqc {
namespace molecule {

/*!
 * \defgroup Molecule Molecule
 *
 * \brief The molecule module contains information about how to make and cluster
 *molecules
 *
 *
 *  The molecule module contains all of the classes which are needed for 
 *  clustering.  Ulitimately everything in the molecule module should support 
 *  the AtomClusterable interface.
 *  
 * @{
 */

/*! \brief Molecule is a class which contains a vector of AtomBasedClusterables
 *
 * At its core molecule is a collection of things which can be collapsed
 * to atoms.  Its main job is allow for clustering of its clusterables.
 *
 */
class Molecule {
  private:
    std::vector<AtomBasedClusterable> elements_;

    Vec3D center_ = {0, 0, 0};
    double mass_ = 0.0;
    int64_t charge_ = 0;

  public:
    Molecule(std::vector<AtomBasedClusterable> c);

    Vec3D const &center() const;
    int64_t charge() const;
    double mass() const;

    std::vector<AtomBasedClusterable>::const_iterator begin() const;
    std::vector<AtomBasedClusterable>::const_iterator end() const;

    int64_t nclusters() const;

    int64_t occupation(unsigned long total_charge) const;

    int64_t core_electrons() const;

    double nuclear_repulsion() const;

};

Molecule read_xyz(std::string const &);

/*!
 * @}
 */

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_MOLECULE_H
