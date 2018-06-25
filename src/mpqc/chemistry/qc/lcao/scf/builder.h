#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_BUILDER_H_

#include <limits>
#include <string>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/expression/formula_registry.h"
#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace lcao {
namespace scf {

/// FockBuilder computes the (2-e part) of the Fock matrix using either the
/// (1-particle) density matrix
/// or the occupied orbital AO coefficients.
template<typename Tile, typename Policy>
class FockBuilder {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  virtual ~FockBuilder() {}

  /// @brief computes the 2-e part of the Fock matrix
  /// @param D the (1-particle) density matrix
  /// @param C the occupied orbital AO coefficient matrix
  /// @param target_precision the target precision
  /// @return the 2-electron contribution to the Fock matrix.
  virtual array_type operator()(
      array_type const &D, array_type const &C,
      double target_precision = std::numeric_limits<double>::epsilon()) = 0;

  virtual void print_iter(std::string const &) = 0;

  virtual void register_fock(const array_type &,
                             FormulaRegistry<array_type> &) = 0;
};

/*!
 * PeriodicFockBuilder computes the (2-e part of) the Fock matrix in periodic HF
 * using the (1-particle) density matrix
 */
template<typename Tile, typename Policy>
class PeriodicFockBuilder {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  virtual ~PeriodicFockBuilder() {}

  /*!
   * \brief This computes the 2-e part of the Fock matrix in periodic HF
   * \param D the (1-particle) density matrix D(μ_0, ν_Rd)
   * \param target_precision the target engine precision
   * \param is_density_diagonal true if D only has a D(μ_0, ν_0) block
   * \return Fock matrix
   */
  virtual array_type operator()(
      array_type const &D,
      double target_precision = std::numeric_limits<double>::epsilon(),
      bool is_density_diagonal = false) = 0;

  virtual void register_fock(const array_type &,
                             FormulaRegistry<array_type> &) = 0;

  /*!
   * \brief This returns the lattice range of Fock representation
   * F(μ_0, ν_R)
   * \return a vector of number of unit cells included in each positive
   * direction
   */
  virtual Vector3i fock_lattice_range() = 0;
};

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_BUILDER_H_
