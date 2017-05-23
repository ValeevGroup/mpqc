#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_PERIODIC_UTIL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_PERIODIC_UTIL_H_

#include "mpqc/chemistry/molecule/unit_cell.h"
#include "mpqc/chemistry/qc/lcao/basis/basis.h"

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor>
    MatrixZ;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> VectorZ;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorD;

// constant
const std::complex<double> I(0.0, 1.0);

namespace mpqc {
namespace lcao {

namespace detail {
/*!
 * \brief This extends 1D tiled range by repeating it multiple times
 *
 * \param tr0 the original TiledRange1 object
 * \param size the number of times for repeating
 * \return the extended TiledRange1 object
 */
TA::TiledRange1 extend_trange1(TA::TiledRange1 const &tr0, int64_t size);

/*!
 * \brief This sorts eigenvalues and eigenvectors
 * in ascending order of the real parts of eigenvalues
 *
 * \param eigVal the vector of complex eigenvalues
 * \param eigVec the complex matrix consisting of complex eigenvectors
 */
void sort_eigen(VectorZ &eigVal, MatrixZ &eigVec);

/*!
 * \brief This takes the ordinal index of a lattice
 * and returns the corresponding direct-space lattice vector
 *
 * \param ord_index the ordinal index of the lattice
 * \param latt_max the range of included lattices
 * \param dcell the direct unit cell params
 * \return the direct-space lattice vector
 */
Vector3d direct_vector(int64_t ord_idx, Vector3i const &latt_max,
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
 * \param input 3D index
 * \param latt_max the range of included lattices
 * \return the ordinal index in direct space
 */
int64_t direct_ord_idx(Vector3i const &in_3D_idx, Vector3i const &latt_max);

/*!
 * \brief This takes the 3D index of a direct lattice
 * and returns the corresponding ordinal index
 *
 * \param x the direct lattice index on x axis
 * \param y the direct lattice index on y axis
 * \param z the direct lattice index on z axis
 * \param latt_max the range of included lattices
 * \return the ordinal index in direct space
 */
int64_t direct_ord_idx(int64_t x, int64_t y, int64_t z,
                       Vector3i const &latt_max);

/*!
 * \brief This takes the ordinal index of a direct lattice
 * and returns the corresponding 3D index
 *
 * \param ord_idx the ordinal index in direct space
 * \param latt_max the range of included lattices
 * \return 3D index in direct space
 */
Vector3i direct_3D_idx(const int64_t ord_idx, Vector3i const &latt_max);

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

namespace gaussian {
namespace detail {

/*!
 * \brief This shifts the origin of a Basis object
 *
 * \note All functions in the basis will be shifted
 * \param basis the original Basis object
 * \param shift the 3D vector of the shift
 * \return the shared pointer of shifted Basis object
 */
std::shared_ptr<Basis> shift_basis_origin(const Basis &basis,
                                          const Vector3d &shift);

/*!
 * \brief This shifts the origin of a Basis object by multiple vectors,
 * and returns a compound Basis that combines all shifted bases
 *
 * \param basis the original Basis object
 * \param shift_base the base position where all shifting vectors start
 * \param nshift the range of included lattices
 * \param dcell the direct unit cell params
 * \return the shared pointer of the compound Basis object
 */
std::shared_ptr<Basis> shift_basis_origin(const Basis &basis,
                                          const Vector3d &shift_base,
                                          const Vector3i &nshift,
                                          const Vector3d &dcell);

}  // namespace detail
}  // namespace gaussian

}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_PERIODIC_UTIL_H_
