#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_LATTICE_UTIL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_LATTICE_UTIL_H_

#include "mpqc/chemistry/molecule/lattice/unit_cell.h"

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace detail {

/*!
 * \brief This takes the ordinal index of a lattice
 * and returns the corresponding direct-space lattice vector
 *
 * \param ord_index the ordinal index of the lattice
 * \param lattice_max the range of included lattices
 * \param dcell the direct unit cell params
 * \return the direct-space lattice vector
 */
Vector3d direct_vector(int64_t ord_idx, Vector3i const &lattice_max,
                       Vector3d const &dcell);

/*!
 * \brief This takes the ordinal index of a k-space (reciprocal-space) lattice
 * and returns the corresponding k-space lattice vector
 *
 * \param ord_idx the ordinal index of the reciprocal lattice
 * \param nk the range of included k points
 * \param dcell the direct unit cell params
 * \return the k-space lattice vector
 */
Vector3d k_vector(int64_t ord_idx, Vector3i const &nk, Vector3d const &dcell);

/*!
 * \brief This takes the 3D index of a direct lattice
 * and returns the corresponding ordinal index
 *
 * \param in_3D_idx input 3D index
 * \param lattice_max the range of included lattices
 * \return the ordinal index in direct space
 */
int64_t direct_ord_idx(Vector3i const &in_3D_idx, Vector3i const &lattice_max);

/*!
 * \brief This takes the 3D index of a direct lattice
 * and returns the corresponding ordinal index
 *
 * \param x the direct lattice index on x axis
 * \param y the direct lattice index on y axis
 * \param z the direct lattice index on z axis
 * \param lattice_max the range of included lattices
 * \return the ordinal index in direct space
 */
int64_t direct_ord_idx(int64_t x, int64_t y, int64_t z,
                       Vector3i const &lattice_max);

/*!
 * \brief This takes the 3D index of a reciprocal lattice
 * and returns the corresponding ordinal index
 *
 * \param in_3D_idx input 3D index
 * \param nk the range of included k points
 * \return the ordinal index in k space
 */
int64_t k_ord_idx(Vector3i const &in_3D_idx, Vector3i const &nk);

/*!
 * \brief This takes the 3D index of a reciprocal lattice
 * and returns the corresponding ordinal index
 *
 * \param x the reciprocal lattice index on x' axis
 * \param y the reciprocal lattice index on y' axis
 * \param z the reciprocal lattice index on z' axis
 * \param nk the range of included k points
 * \return the ordinal index in k space
 */
int64_t k_ord_idx(int64_t x, int64_t y, int64_t z, Vector3i const &nk);

/*!
 * \brief This takes the ordinal index of a direct lattice
 * and returns the corresponding 3D index
 *
 * \param ord_idx the ordinal index in direct space
 * \param lattice_max the range of included lattices
 * \return 3D index in direct space
 */
Vector3i direct_3D_idx(const int64_t ord_idx, Vector3i const &lattice_max);

/*!
 * \brief This shifts the position of a Molecule object
 *
 * \note All atom positions will be shifted
 * \param mol the Molecule object
 * \param shift the 3D vector of the shift
 * \return the shared pointer of the shifted Molecule object
 */
std::shared_ptr<Molecule> shift_mol_origin(Molecule const &mol,
                                           Vector3d const &shift);

}  // namespace detail
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_LATTICE_UTIL_H_
