#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_UNITCELL_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_UNITCELL_H_

#include "./molecule.h"

#include <iosfwd>
#include <vector>

#include <Eigen/Dense>

#include "mpqc/util/keyval/keyval.hpp"

namespace mpqc {
class UnitCell : public Molecule {
 private:
  Vector3i R_max_ = {
      0, 0, 0};  // range of expansion of Bloch Gaussians in AO Gaussians
  Vector3i RJ_max_ = {0, 0, 0};       // range of Coulomb operation
  Vector3i RD_max_ = {0, 0, 0};       // range of density representation
  Vector3i nk_ = {1, 1, 1};           // # of k points in each direction
  Vector3d dcell_ = {0.0, 0.0, 0.0};  // direct unit cell params (in a.u.)
 public:
  UnitCell() = default;

  /** \brief KeyVal constructor for Periodic System
   *
   *  <table border="1">
   *
   *  <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>
   *
   *  <tr><td><tt>file_name</tt><td>string<td>none<td>This gives
   *    the name of a XYZ file, from which the nuclear coordinates will be
   *    read (the XYZ format is described
   *    <a href="http://en.wikipedia.org/wiki/XYZ_file_format">here</a>).
   *
   *  <tr><td><tt>charge</tt><td>int<td>0<td> the charge of this molecule
   *
   *  <tr><td><tt>sort_input</tt><td>boolean<td>true<td>If true, atoms
   *    will be resorted based on their distance from the center of mass.
   *
   *  <tr><td><tt>sort_origin</tt><td>boolean<td>false<td>If true, sort atoms
   *  from origin {0.0, 0.0, 0.0}
   *
   *  <tr><td><tt>n_cluster</tt><td>int<td>0<td> If nonzero, cluster moleucle by
   * n_cluster
   *
   *  <tr><td><tt>attach_hydrogen</tt><td>bool<td>true<td> use
   * attach_hydrogen_kmeans when clustering
   *
   *  <tr><td><tt>direct_lattice_vector</tt><td>int<td>0<td> This gives three
   * direct lattice vectors:
   *  [Vx, 0.0, 0.0], [0.0, Vy, 0.0], and [0.0, 0.0, Vz]. Units are Bohr.
   *
   *  <tr><td><tt>max_lattice_sum</tt><td>int<td>0<td> Go out [Nx, Ny, Nz] cells
   *  in both positive and negative directions of three dimensions. Applied in
   * all lattice sums
   *
   *  <tr><td><tt>k_points</tt><td>int<td>0<td> number of k points
   *  </table>
   *
   *  example input:
   *  \code
   *  "molecule": {
   *    "type": "UnitCell",
   *    "charge": 0,
   *    "file_name": "water.xyz",
   *    "sort_input": true,
   *    "n_cluster": 20
   *    "max_lattice_sum": [0, 0, 10],
   *    "direct_lattice_vector": [0.0, 0.0, 2.8],
   *    "k_points": [1, 1, 10]
   *  }
   *  \endcode
   *
   */
  UnitCell(const KeyVal& kv);

  /// Return the nuclear repulsion energy of the Periodic System.
  double nuclear_repulsion(Vector3i RJ_max) const;

  ~UnitCell() = default;

  /// Print out molecule and atom information
  void print(std::ostream& out) const override;

  /// Return dcell_
  Vector3d dcell() { return dcell_; }
  Vector3i R_max() { return R_max_; }
  Vector3i RD_max() { return RD_max_; }
  Vector3i RJ_max() { return RJ_max_; }
  Vector3i nk() { return nk_; }
};

}  // mpqc namespace

#endif  // SRC_MPQC_CHEMISTRY_MOLECULE_UNITCELL_H
