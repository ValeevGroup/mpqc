#pragma once
#ifndef MPQC_MOLECULE_MOLECULE_H
#define MPQC_MOLECULE_MOLECULE_H

#include <mpqc/chemistry/molecule/molecule_fwd.h>
#include <mpqc/chemistry/molecule/atom_based_cluster.h>
#include <mpqc/util/keyval/keyval.hpp>

#include <iosfwd>
#include <vector>


namespace mpqc {
namespace molecule {

/*!
 * \defgroup Molecule Molecule
 *
 * \brief The molecule module contains information about how to make and cluster
 *molecules
 *
 * @{
 */

/*! \brief Molecule is a class which contains a vector of AtomBasedClusterables
 *
 * At its core molecule is a collection of clusters that can be collapsed
 * to atoms.  Its main job is allow for clustering of its clusterables.
 *
 */
class Molecule : public DescribedClass {
  private:
    std::vector<AtomBasedClusterable> elements_;

    Vec3D com_ = {0, 0, 0}; /// Center of Mass
    double mass_ = 0.0;
    int64_t charge_ = 0; /// Net charge (# protons - # electrons)

    void init(std::istream &file, bool sort_input);

  public:
    Molecule() = default;

    /*! \brief KeyVal constructor for Molecule
     *
     *  The KeyVal constructor expects a field called file_name and an optional
     *  input sort_input, which are the name of a file that contains molecular
     *  coordinates and a boolean telling the constructor whether the atoms
     *  should be sorted based on their distance from the center of mass. By
     *  default sort_input is set to true.
     *
     *  The clustering of the KeyVal constructor is currently one cluster per
     *  atom, but this may change as KeyVal input is finalized.
     */
    Molecule(const KeyVal &kv);

    /*! \brief Constructor to build Molecule from stream.
     *
     * This constructor has same parameters as the KeyVal constructor.
     */
    Molecule(std::istream &file_stream, bool sort_input = true);

    /*! \brief Constructor to build Molecule from a vector of clusterables
     *
     * This constructor takes a vector of AtomBasedClusterables and uses it to
     * initialize the Molecule's Clusterables. If sort_input is true the
     * Clusterables are sorted based on distance from the center of mass.
     */
    Molecule(std::vector<AtomBasedClusterable> c, bool sort_input = true);

    /// Function to set the mass of the Molecule
    void set_mass(double mass) { Molecule::mass_ = mass; }

    /// Function to set the charge of the Molecule
    void set_charge(int64_t charge) { Molecule::charge_ = charge; }

    /// Charge of the Molecule: charge = (# protons - # electrons)
    int64_t charge() const { return charge_; }

    /// Mass of the Molecule
    double mass() const { return mass_; }

    /// Iterator to the first clusterable in the Molecule
    std::vector<AtomBasedClusterable>::const_iterator begin() const {
        return elements_.begin();
    }

    /// Iterator to one past the final clusterable in the Molecule
    std::vector<AtomBasedClusterable>::const_iterator end() const {
        return elements_.end();
    }

    /// Number of clusters in the Molecule
    int64_t nclusters() const { return elements_.size(); }

    /// Vector containing the Clusterable that make up the Molecule
    std::vector<AtomBasedClusterable> const &clusterables() const {
        return elements_;
    }

    /*! \brief Returns the difference between the Molecule's charge and the
     * charge passed into the function.
     *
     * \note this interface is a bit weird and will likely change with updates
     * to KeyVal construction.
     */
    int64_t occupation(unsigned long total_charge) const {
        return charge_ - total_charge;
    }

    /// Computes the number of core electrons in the Molecule.
    int64_t core_electrons() const;

    /// Returns the nuclear repulsion energy of the Molecule.
    double nuclear_repulsion() const;

    /*! \brief A vector of all atoms in the Molecule
     *
     * This function will loop over the clusterables and extract their atoms in
     * a recursive fashion.
     *
     * \note Returns copies of the atoms, not a reference to the atoms stored
     * in the Clusterables of the Molecule.
     */
    std::vector<Atom> atoms() const;

    /*! \brief Center of mass of the Molecule.
     *
     * Necessary to satisfy the AtomBasedClusterable interface so Molecules are
     * also clusterable.
    */
    Vec3D const &com() const { return com_; }
};

/// Make Molecules printable
std::ostream &operator<<(std::ostream &, Molecule const &);

/*! @} */

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_MOLECULE_H
