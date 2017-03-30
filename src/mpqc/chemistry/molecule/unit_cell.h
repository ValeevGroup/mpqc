#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_UNIT_CELL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_UNIT_CELL_H_

#include "mpqc/chemistry/molecule/molecule.h"

#include <iosfwd>
#include <vector>

#include <Eigen/Dense>

#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
class UnitCell : public Molecule {
 private:
  /// \todo extend to nonorthogonal unit cells
  Vector3d dcell_ = {0.0, 0.0, 0.0};  ///< direct unit cell params (in a.u.)

 public:
  UnitCell() = default;

  /** \brief KeyVal constructor for Periodic System
   *
   *  <table border="1">
   *
   *  <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>
   *
   *  <tr><td><tt>type</tt><td>int<td>0<td> the type of this molecule. If type
   *    is UnitCell, periodic calculations can be performed
   *
   *  <tr><td><tt>file_name</tt><td>string<td>none<td>This gives
   *    the name of a XYZ file, from which the nuclear coordinates will be
   *    read (the XYZ format is described
   *    <a href="http://en.wikipedia.org/wiki/XYZ_file_format">here</a>).
   *
   *  <tr><td><tt>charge</tt><td>int<td>0<td> the charge of unitcell
   *
   *  <tr><td><tt>sort_input</tt><td>boolean<td>true<td>If true, atoms
   *    will be resorted based on their distance from the center of mass.
   *
   *  <tr><td><tt>sort_origin</tt><td>boolean<td>false<td>If true, sort atoms
   *    from origin {0.0, 0.0, 0.0}
   *
   *  <tr><td><tt>n_cluster</tt><td>int<td>0<td> If nonzero, cluster moleucle by
   *    n_cluster
   *
   *  <tr><td><tt>attach_hydrogen</tt><td>bool<td>true<td> use
   *    attach_hydrogen_kmeans when clustering
   *
   *  <tr><td><tt>lattice_param</tt><td>int<td>0<td> This gives lattice
   *    parameters (orthohombic, tetragonal, cubic lattice only)
   *
   *  </table>
   *
   *  example input:
   *  \code
   *  "unitcell": {
   *    "type": "UnitCell",
   *    "charge": 0,
   *    "file_name": "water.xyz",
   *    "sort_input": true,
   *    "lattice_param": [0.0, 0.0, 2.672359],
   *  }
   *  \endcode
   *
   */
  UnitCell(const KeyVal &kv);

  /*!
   * \brief This computes the nuclear repulsion energy of the unit cell with other cell within a range.
   * \note this includes the intra-cell repulsion.
   * \param RJ_max the range of nuclear repulsion interaction; all cells with [- \c RJ_max .. RJ_max] are included
   * \return nuclear repulsion energy
   */
  double nuclear_repulsion_energy(Vector3i RJ_max) const;

  ~UnitCell() { }

  /*!
   * \brief This returns direct unit cell params
   * \return dcell_
   */
  Vector3d dcell() const { return dcell_; }
};

/// Make UnitCell printable
std::ostream &operator<<(std::ostream &, UnitCell const &);

}  // mpqc namespace

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_UNIT_CELL_H_
