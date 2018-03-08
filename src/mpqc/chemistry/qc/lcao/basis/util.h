#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_SHIFT_BASIS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_SHIFT_BASIS_H_

#include "mpqc/chemistry/qc/lcao/basis/basis.h"

#include "mpqc/chemistry/molecule/lattice/util.h"

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

using func_offset_list = std::unordered_map<size_t, std::tuple<size_t, size_t>>;

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
 * \param is_half_range true if only half range of \c nshift is needed, false
 * if the full range of \c nshift is needed.
 * \return the shared pointer of the compound Basis object
 */
std::shared_ptr<Basis> shift_basis_origin(const Basis &basis,
                                          const Vector3d &shift_base,
                                          const Vector3i &nshift,
                                          const Vector3d &dcell,
                                          const bool is_half_range = false);

/*!
 * \brief This computes shell offsets for every cluster in a basis
 * \param basis
 * \return a list of <key, mapped value> pairs with
 * key: cluster index
 * mapped value: index of first shell in a cluster
 */
std::unordered_map<size_t, size_t> compute_shell_offset(const Basis &basis);

/*!
 * \brief This computes cluster & basis function offsets for every shell in a
 * cluster.
 * \param cluster a cluster (a.k.a. std::vector<Shell>)
 * \param bf_first basis function index of the first function in this \c
 * cluster
 *
 * \return a list of <key, mapped value> pairs with
 * key: shell index
 * mapped value: {cluster function offset, basis function offset} tuple
 */
func_offset_list compute_func_offset_list(const ShellVec &cluster,
                                          const size_t bf_first);

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_SHIFT_BASIS_H_
