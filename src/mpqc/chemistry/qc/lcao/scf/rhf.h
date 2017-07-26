//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RHF_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/expression/orbital_space.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/density_builder.h"
#include "mpqc/chemistry/qc/lcao/wfn/ao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

/// Closed-Shell (Spin-)Restricted SCF
template <typename Tile, typename Policy>
class RHF
    : public AOWavefunction<Tile, Policy>,
      public Provides<Energy,
                      CanonicalOrbitalSpace<TA::DistArray<Tile, Policy>>,
                      PopulatedOrbitalSpace<TA::DistArray<Tile, Policy>>> {
 public:
  using array_type = TA::DistArray<Tile, Policy>;

  RHF() = default;

  // clang-format off
  /**
   * KeyVal constructor for RHF
   *
   * keywords: takes all keywords from AOWavefunction
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | max_iter | int | 30 | maximum number of iteration |
   * | density_builder | string | eigen_solve | type of DensityBuilder, valid values are \c eigen_solve (use ESolveDensityBuilder) and \c purification (use PurificationDensityBuilder) |
   * | localize | bool | false | if localize in DensityBuilder |
   * | t_cut_c | double | 0.0 | threshold in DensityBuilder, SparsePolicy only |
   * | decompo_type | string | conditioned | (cholesky_inverse, inverse_sqrt, conditioned) only valid if use ESolveDensityBuilder |
   * | s_tolerance | double | 1.0e8 | S condition number threshold in DensityBuilder, valid when decompo_type is set to conditioned |
   *
   */
  // clang-format on
  RHF(const KeyVal& kv);

  virtual ~RHF() { }

  void obsolete() override;

  /// @return the number of electrons (twice the number of occupied orbitals,
  /// hence an even integer)
  size_t nelectrons() const { return nelectrons_; }

  // this is almost evil
  inline void set_fock(array_type f) { F_ = f; }

  /// these implement Energy::Provider methods
  bool can_evaluate(Energy* result) override;
  void evaluate(Energy* result) override;

  // these implement CanonicalOrbitalSpace::Provider methods

  /// @return true if using eigen solver and orbitals are not localized
  bool can_evaluate(CanonicalOrbitalSpace<array_type>* = nullptr) override;
  /// Computes all canonical orbitals, annotated with orbital energies
  void evaluate(CanonicalOrbitalSpace<array_type>* result,
                double target_precision, std::size_t target_blocksize) override;

  // these implement PopulatedOrbitalSpace::Provider methods

  /// @return true if using eigen solver
  bool can_evaluate(PopulatedOrbitalSpace<array_type>* = nullptr) override;
  /// Computes all or occupied orbitals, annotated with occupancies
  void evaluate(PopulatedOrbitalSpace<array_type>* result,
                double target_precision, std::size_t target_blocksize) override;

 protected:
  double energy_;
  std::size_t max_iter_;

  std::size_t nelectrons_;

  array_type H_;
  array_type S_;
  array_type F_;
  array_type F_diis_;
  array_type D_;
  array_type C_;  //!< occupied orbitals only

  std::unique_ptr<scf::FockBuilder<Tile, Policy>> f_builder_;
  std::unique_ptr<scf::DensityBuilder<Tile, Policy>> d_builder_;

  std::vector<double> rhf_times_;
  std::vector<double> d_times_;
  std::vector<double> build_times_;

  std::string density_builder_str_;
  bool localize_;
  double t_cut_c_;

 private:
  // to expose these need to wrap into if_computed
  inline array_type const& overlap() const { return S_; }
  inline array_type const& fock() const { return F_; }
  inline array_type const& density() const { return D_; }
  inline array_type const& coefficents() const { return C_; }
  inline const double energy() const { return energy_; }

  double compute_energy() const;

  /** Function to compute the density to the desired accuracy.
   *
   * Takes some form of integral and does the closed-shell scf iterations.  The
   * place to
   * specialize is in build_fock.
   *
   * @param max_iters the maximum number of iterations
   * @param target_precision the target SCF convergence threshold; SCF is
   * converged if
   *        the relative energy change, \f$ |E_i - E_{i-1}|/E_i \f$, and the
   *        per-element AO-basis orbital gradient, \f$
   * ||\mathbf{F},\mathbf{P}||_2 / n^2 \f$,
   *        are below \c thresh .
   * @throws MaxIterExceeded if exceeded the limit on the number of iteration
   */
  void solve(int64_t max_iters, double target_precision);
  double computed_precision_ = std::numeric_limits<double>::max();

  /// does compute-heavy initialization (integrals, guess) and runs SCF
  /// @param target_precision the target SCF convergence threshold (see
  /// RHF::solve() )
  void do_evaluate(double target_precision);

  void init(const KeyVal& kv);
  virtual void init_fock_builder();
  void compute_density();
  void build_F();

  const KeyVal kv_;
};

/**
 *
 * RIRHF class, fock_builder uses density fitting
 */
template <typename Tile, typename Policy>
class RIRHF : public RHF<Tile, Policy> {
 public:
  RIRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

/**
 * CadfRHF class, using direct traditional density fitting for J and the
 * Concentric Atomic Density Fitting Approach for K.
 *
 */
template <typename Tile, typename Policy>
class CadfRHF : public RHF<Tile, Policy> {
 public:
  /*!
   * Parameter tcutc can be set to truncate elements of the molecular orbitals,
   * by default it is 0.0 ensuring no truncation
   *
   * A further approximation called force shape may be applied by including the
   * force_shape_threshold keyword and assigning a value greater than 0 to it.
   *
   * Finally if force_shape_threshold != 0 then tcutc will be defaulted to 1e-4,
   * but will still be settable by the user.
   */
  CadfRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;

  double force_shape_threshold_ = 0.0;
  double tcutc_ = 0.0;
  bool secadf_ = false;
  bool aaab_ = false;
};

/**
 * DirectRIRHF, fock_builder uses density fitting with 3-center integrals
 * computed on the fly
 */
template <typename Tile, typename Policy>
class DirectRIRHF : public RHF<Tile, Policy> {
 public:
  DirectRIRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

/**
 * DirectRHF, fock_builder uses 4-center integrals computed on the fly.
 */
template <typename Tile, typename Policy>
class DirectRHF : public RHF<Tile, Policy> {
 public:
  DirectRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

/**
 * RIJ-EXACTK-RHF, fock_builder is overide to use direct four center integral
 */
template <typename Tile, typename Policy>
class RIJEXACTKRHF : public RHF<Tile, Policy> {
 public:
  RIJEXACTKRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

#if TA_DEFAULT_POLICY == 0
extern template class RHF<TA::TensorD, TA::DensePolicy>;
extern template class RIRHF<TA::TensorD, TA::DensePolicy>;
extern template class DirectRHF<TA::TensorD, TA::DensePolicy>;
extern template class DirectRIRHF<TA::TensorD, TA::DensePolicy>;
extern template class CadfRHF<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class RHF<TA::TensorD, TA::SparsePolicy>;
extern template class RIRHF<TA::TensorD, TA::SparsePolicy>;
extern template class DirectRHF<TA::TensorD, TA::SparsePolicy>;
extern template class DirectRIRHF<TA::TensorD, TA::SparsePolicy>;
extern template class CadfRHF<TA::TensorD, TA::SparsePolicy>;
extern template class RIJEXACTKRHF<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "rhf_impl.h"
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RHF_H_
