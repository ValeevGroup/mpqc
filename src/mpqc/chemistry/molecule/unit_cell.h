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

  // clang-format off
  /** \brief KeyVal constructor for a unit cell of an orthohombic, tetragonal, or cubic lattice.
   *
   * \param kv The KeyVal object will be queried for all keywords of OrbitalBasisRegistry, and the following keywords:
   *  | Keyword | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |\c lattice_param | int | 0 | This gives unit cell dimensions. |
   *
   *  example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "unitcell": {
   *    "file_name": "water.xyz",
   *    "sort_input": true,
   *    "lattice_param": [0.0, 0.0, 2.672359],
   *  }
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  // clang-format on
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
