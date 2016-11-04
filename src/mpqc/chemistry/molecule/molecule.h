#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_MOLECULE_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_MOLECULE_H_

#include <iosfwd>
#include <vector>

#include "./atom_based_cluster.h"
#include "./molecule_fwd.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {

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
 protected:
  std::vector<AtomBasedClusterable> elements_;

  Vector3d com_ = {0, 0, 0};  /// Center of Mass
  double mass_ = 0.0;
  int64_t charge_ = 0;        /// Net charge (# protons - # electrons)
  int64_t total_charge_ = 0;  // total charge # protons

  void init(std::istream &file, bool sort_input);

  void init(std::istream &file, Vector3d const &point);

 public:
  Molecule() = default;

  /** \brief KeyVal constructor for Molecule
   *
   *
   *  | KeyWord | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |file_name|string|none|This gives the name of a XYZ file, from which the
   * nuclear coordinates will be read |
   *  |||| (the XYZ format is described <a
   * href="http://en.wikipedia.org/wiki/XYZ_file_format">here</a>).|
   *  |charge|int|0|the charge of this molecule|
   *  |sort_input|boolean|true|If true, sort atoms from origin {0.0, 0.0, 0.0} |
   *  |sort_origin|boolean|false|sort atoms from origin {0.0, 0.0, 0.0} |
   *  |n_cluster|int|0|If nonzero, cluster moleucle by n_cluster|
   *  |attach_hydrogen|boolean|true|use attach_hydrogen_kmeans when clustering|
   *
   *
   *
   *  example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "molecule": {
   *    "type": "Molecule",
   *    "charge": 0,
   *    "file_name": "water.xyz",
   *    "sort_input": true,
   *    "n_cluster": 20
   *  }
   *
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  Molecule(const KeyVal &kv);

  /*! \brief Constructor to build Molecule from stream.
   *
   * This constructor has same parameters as the KeyVal constructor.
   */
  Molecule(std::istream &file_stream, bool sort_input = true);

  /*! \brief Constructor to build Molecule from stream.
   *
   * This constructor changes the point from which the molecule is sorted
   */
  Molecule(std::istream &file_stream, Vector3d const &point);

  /*! \brief Constructor to build Molecule from a vector of clusterables
   *
   * This constructor takes a vector of AtomBasedClusterables and uses it to
   * initialize the Molecule's Clusterables. If sort_input is true the
   * Clusterables are sorted based on distance from the center of mass.
   */
  Molecule(std::vector<AtomBasedClusterable> c, bool sort_input = true);

  ~Molecule();

  /*! \brief A function to sort the molcule's clusters from a given point.
   */
  void sort_from_point(Vector3d const &point);

  /// Function to set the mass of the Molecule
  void set_mass(double mass) { Molecule::mass_ = mass; }

  /// Function to set the charge of the Molecule
  void set_charge(int64_t charge) { Molecule::charge_ = charge; }

  /// Charge of the Molecule: charge = (# protons - # electrons)
  int64_t charge() const { return charge_; }

  /// total nuclear charge # of protons
  int64_t total_charge() const { return total_charge_; }

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

  /*! \brief Returns the difference between the Molecule's total_charge and the
   * charge passed into the function.
   */
  int64_t occupation(unsigned long charge) const {
    return total_charge_ - charge;
  }

  /// return the occupation number
  int64_t occupation() const { return total_charge_ - charge_; }

  /// Computes the number of core electrons in the Molecule.
  int64_t core_electrons() const;

  /// Returns the nuclear repulsion energy of the Molecule.
  virtual double nuclear_repulsion() const;

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
  Vector3d const &com() const { return com_; }

  /// Print Molecule information
  virtual void print(std::ostream & os) const;
};

/// Make Molecules printable
std::ostream &operator<<(std::ostream &, Molecule const &);

/*! @} */

}  // namespace mpqc

#endif  // MPQC_MOLECULE_MOLECULE_H
