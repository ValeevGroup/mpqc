//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RHF_H_

#include <tiledarray.h>

#include "mpqc/mpqc_config.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/chemistry/qc/scf/builder.h"
#include "mpqc/chemistry/qc/scf/density_builder.h"
#include "mpqc/chemistry/qc/wfn/ao_wfn.h"
/**
 *  RHF Class of AOWfn
 *
 */
namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class RHF : public AOWavefunction<Tile, Policy>, public CanEvaluate<Energy> {
 public:
  using array_type = TA::DistArray<Tile,Policy>;

  RHF() = default;

  /**
   * KeyVal constructor for RHF
   *
   * keywords: takes all keywords from AOWavefunction
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | max_iter | int | 30 | maximum number of iteration |
   * | density_builder | string | eigen_solve | type of DensityBuilder (eigen_solve->ESolveDensityBuilder, purification->PurificationDensityBuilder) |
   * | localize | bool | false | if localize in DensityBuilder |
   * | t_cut_c | double | 0.0 | threshold in DensityBuilder, SparsePolicy only |
   * | decompo_type | string | cholesky_inverse | (cholesky_inverse, inverse_sqrt, conditioned_inverse) only valid if use ESolveDensityBuilder |
   *
   */

  RHF(const KeyVal& kv);

  virtual ~RHF() = default;

  void obsolete() override;

  inline void set_fock(array_type f) {F_ = f;}
 protected:
  double energy_;
  std::size_t max_iter_;

  array_type H_;
  array_type S_;
  array_type F_;
  array_type F_diis_;
  array_type D_;
  array_type C_;

  std::unique_ptr<scf::FockBuilder<Tile,Policy>> f_builder_;
  std::unique_ptr<scf::DensityBuilder<Tile,Policy>> d_builder_;

  std::vector<double> rhf_times_;
  std::vector<double> d_times_;
  std::vector<double> build_times_;

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
   * Takes some form of integral and does the closed-shell scf iterations.  The place to
   * specialize is in build_fock.
   *
   * @param max_iters the maximum number of iterations
   * @param thresh the target SCF convergence threshold; SCF is converged if
   *        the relative energy change, \f$ |E_i - E_{i-1}|/E_i \f$, and the
   *        per-element AO-basis orbital gradient, \f$ ||\mathbf{F},\mathbf{P}||_2 / n^2 \f$,
   *        are below \c thresh .
   * @throws MaxIterExceeded if exceeded the limit on the number of iteration
   */
  void solve(int64_t max_iters, double thresh);

  void init(const KeyVal& kv);
  virtual void init_fock_builder();
  void compute_density();
  void build_F();

  bool can_evaluate(Energy* energy) override;
  void evaluate(Energy* result) override;

  const KeyVal kv_;
};

/**
 *
 * RIRHF class, fock_builder is overide to use three center integral
 */
template <typename Tile, typename Policy>
class RIRHF : public RHF<Tile,Policy> {
 public:
  RIRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

/**
 * DirectRIRHF, fock_builder is overide to use direct three center integral
 */
template <typename Tile, typename Policy>
class DirectRIRHF : public RHF<Tile,Policy> {
 public:
  DirectRIRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

/**
 * DirectRHF, fock_builder is overide to use direct four center integral
 */
template <typename Tile, typename Policy>
class DirectRHF : public RHF<Tile,Policy> {
 public:
  DirectRHF(const KeyVal& kv);

 private:
  void init_fock_builder() override;
};

#if TA_DEFAULT_POLICY == 0
extern template class RHF<TA::TensorD, TA::DensePolicy>;
extern template class RIRHF<TA::TensorD, TA::DensePolicy>;
extern template class DirectRHF<TA::TensorD, TA::DensePolicy>;
extern template class DirectRIRHF<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class RHF<TA::TensorD, TA::SparsePolicy>;
extern template class RIRHF<TA::TensorD, TA::SparsePolicy>;
extern template class DirectRHF<TA::TensorD, TA::SparsePolicy>;
extern template class DirectRIRHF<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#include "rhf_impl.h"
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RHF_H_
