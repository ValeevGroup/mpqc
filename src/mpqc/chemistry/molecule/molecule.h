#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_MOLECULE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_MOLECULE_H_

#include <iosfwd>
#include <vector>

#include "./atom_based_cluster.h"
#include "./molecule_fwd.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/observer.h"

namespace mpqc {

/*!
 * \addtogroup ChemistryMolecule
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
class Molecule : virtual public DescribedClass,
                 public utility::Observable<Molecule> {
 private:
  std::vector<AtomBasedClusterable> elements_;

  Vector3d com_ = {0, 0, 0};  /// Center of Mass
  double mass_ = 0.0;
  int64_t total_atomic_number_ = 0;  // sum of atomic numbers of all atoms
  double natoms_ = 0.0;

  /// Initializes the molecule using a vector of atoms.
  /// @note this does not recluster.
  /// @param atoms an rvalue ref to a vector of Atom objects opaqued as AtomBasedClusterable objects
  /// @param sort_input if true, will resort atoms in the order of increasing distance from
  ///        a reference point.
  /// @param ref_point the reference point to use for sorting. By default will use the center
  ///        of mass.
  void init_atoms(std::vector<AtomBasedClusterable> &&atoms, bool sort_input,
            const Vector3d *const ref_point = nullptr);

  /// @brief Reads atoms from an XYZ file.
  /// @return the vector of Atom objects opaqued as AtomBasedClusterable objects
  static std::vector<AtomBasedClusterable> read_xyz(std::istream &file);

 public:
  /// the only way to mutate coordinates is via MolecularCoordinates
  friend class MolecularCoordinates;
  void update(const std::vector<Atom> &atoms);

 public:
  Molecule() = default;

  // clang-format off
  /** \brief The KeyVal constructor.
   *  \param kv the KeyVal object. The following keywords will be queried in \c kv :
   *
   *  | Keyword | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |\c file_name|string|none|This gives the name of a XYZ file, from which the nuclear coordinates will be read (the XYZ format is described <a href="http://en.wikipedia.org/wiki/XYZ_file_format">here</a>). Only rank 0 in the current madness::World will read the file.|
   *  |\c atoms|array|none|Will query this if \c file_name not given. Each element of the array must specify an Atom object (see the KeyVal ctor of Atom for more info).|
   *  |\c units|string|angstrom|If \c atoms given, will use this to determine the units in which the coordinates are given. Valid choices are \c "angstrom" and \c "bohr" See UnitFactory for more information.|
   *  |\c sort_input|boolean|true|If true, sort atoms in the order of increasing distance from the center of mass |
   *  |\c n_cluster|int|0|If nonzero, divide the Molecule into \c n_cluster clusters|
   *  |\c attach_hydrogen|boolean|true|use attach_hydrogen_kmeans when clustering|
   *
   *  example input 1:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "molecule": {
   *    "file_name": "water20.xyz",
   *    "n_cluster": 20
   *  }
   *
   * ~~~~~~~~~~~~~~~~~~~~~
   *
   *  example input 2:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "molecule": {
   *    "atoms" : [
   *      { "element" : "O", "xyz" : [0.0,  0.0, 0.0] },
   *      { "element" : "H", "xyz" : [0.0,  1.0, 1.0] },
   *      { "element" : "H", "xyz" : [0.0, -1.0, 1.0] }
   *    ]
   *    "n_cluster": 1
   *  }
   *
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  // clang-format on
  Molecule(const KeyVal &kv);

  /*! \brief Constructor to build Molecule from an XYZ-format stream.
   * \note This is a non-clustering ctor, i.e. produces a Molecule consisting of 1 cluster.
   */
  Molecule(std::istream &file_stream, bool sort_input = true);

  /*! \brief Constructor to build Molecule from stream.
   * \note This is a non-clustering ctor, i.e. produces a Molecule consisting of 1 cluster.
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

  /// @return the sum of atomic numbers of all atoms
  /// @sa Atom::atomic_number()
  int64_t total_atomic_number() const { return total_atomic_number_; }

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

  /// Computes the number of core electrons in the Molecule.
  int64_t core_electrons() const;

  /// @return the nuclear repulsion energy of the Molecule.
  double nuclear_repulsion_energy() const;

  /*! \brief A vector of all atoms in the Molecule
   *
   * This function will loop over the clusterables and extract their atoms in
   * a recursive fashion.
   *
   * \note Returns copies of the atoms, not a reference to the atoms stored
   * in the Clusterables of the Molecule.
   */
  std::vector<Atom> atoms() const;

  /// @return the number of atoms returned by atoms()
  size_t natoms() const;

  /*! \brief Center of mass of the Molecule.
   *
   * Necessary to satisfy the AtomBasedClusterable interface so Molecules are
   * also clusterable.
  */
  Vector3d const &com() const { return com_; }
};

/// Make Molecules printable
std::ostream &operator<<(std::ostream &, Molecule const &);

/*! @} */

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_MOLECULE_H_
