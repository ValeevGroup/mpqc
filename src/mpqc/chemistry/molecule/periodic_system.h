#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_PERIODIC_SYSTEM_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_PERIODIC_SYSTEM_H_

#include <iosfwd>
#include <vector>

#include "./molecule.h"
#include <mpqc/util/keyval/keyval.hpp>

#include <Eigen/Dense>

typedef Eigen::Vector3i Vec3I;
typedef Eigen::Vector3d Vec3D;

namespace mpqc {
namespace molecule {

class PeriodicSystem : public Molecule {
 private:
  Vec3I R_max_ = {0, 0, 0};  // range of expansion of Bloch Gaussians in AO Gaussians
  Vec3I RJ_max_ = {0, 0, 0};       // range of Coulomb operation
  Vec3I RD_max_ = {0, 0, 0};       // range of density representation
  Vec3I nk_ = {1, 1, 1};           // # of k points in each direction
  Vec3D dcell_ = {0.0, 0.0, 0.0};  // direct unit cell params (in a.u.)
 public:
  PeriodicSystem() = default;

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
   *  </table>
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
   *
   *  example input:
   *  "molecule": {
   *    "type": "PeriodicSystem",
   *    "charge": 0,
   *    "file_name": "water.xyz",
   *    "sort_input": true,
   *    "n_cluster": 20
   *    "max_lattice_sum": [0, 0, 10],
   *    "direct_lattice_vector": [0.0, 0.0, 2.8],
   *    "k_points": [1, 1, 10]
   *  }
   *
   */
  PeriodicSystem(const KeyVal& kv);

  /// Return the nuclear repulsion energy of the Periodic System.
  double nuclear_repulsion() const override;

  ~PeriodicSystem() = default;

  /// Print out molecule and atom information
  void print(std::ostream& out) const override;

  /// Return dcell_
  Vec3D dcell() {return dcell_;}
  Vec3I R_max() {return R_max_;}
  Vec3I RD_max() {return RD_max_;}
  Vec3I RJ_max() {return RJ_max_;}
  Vec3I nk() {return nk_;}

};

}  // molecule namespace
}  // mpqc namespace

#endif  // SRC_MPQC_CHEMISTRY_MOLECULE_PERIODIC_SYSTEM_H
