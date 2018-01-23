#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"

#include <memory>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/wfn/ao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"

// constant: imaginary unit i
const std::complex<double> I(0.0, 1.0);

namespace mpqc {
namespace lcao {

using MatrixzVec = std::vector<MatrixZ>;
using VectorzVec = std::vector<VectorZ>;
using VectordVec = std::vector<VectorD>;
using Matrix = RowMatrixXd;

/**
 * complex-valued Restricted Hartree-Fock class
 */
template <typename Tile, typename Policy>
class zRHF : public PeriodicAOWavefunction<Tile, Policy>,
             public Provides<Energy /*,
          CanonicalOrbitalSpace<TA::DistArray<TA::TensorZ, TA::SparsePolicy>>,
          PopulatedOrbitalSpace<TA::DistArray<TA::TensorZ, TA::SparsePolicy>>*/> {
 public:
  using array_type = typename PeriodicAOWavefunction<Tile, Policy>::ArrayType;
  using factory_type =
      typename PeriodicAOWavefunction<Tile, Policy>::AOIntegral;
  using array_type_z = TA::DistArray<TA::TensorZ, Policy>;

  zRHF() = default;

  // clang-format off
  /**
   * \brief KeyVal constructor for zRHF
   *
   * \param kv The KeyVal object will take all keywords from PeriodicAOWavefunction, and the following keywords
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | max_iter | int | 30 | maximum number of iteration |
   * | soad_guess | bool | true | if use SOAD guess for initial Fock build |
   * | print_detail | bool | false | if print extra computation&time info |
   * | max_condition_num | double | 1.0e8 | maximum condition number for overlap matrix |
   * | k_points | array<int, 3> | none | number of k points in each direction of the first Brillouin zone |
   * | print_max_item | int | 100 | maximum number of items/lines that can be printed in the list of condition numbers |
   * | fock_mixing | double | 0.0 | mixing of Fock matrices in reciprocal space |
   * | diis | string | none | the choice of DIIS method: none, gamma_point, all_k, sloshing |
   * | diis_start | unsigned int | 1 | the DIIS extrapolation will begin on the iteration given by this integer |
   * | diis_num_vecs | unsigned int | 5 | maximum number of data sets to store |
   * | diis_damping | double | 0.0 | this nonnegative floating point number is used to dampen the DIIS extrapolation |
   * | diis_mixing | double | 0.0 | this nonnegative floating point number is used to dampen the DIIS extrapolation by mixing the input Fock with the output Fock for each iteration |
   * | diis_num_iters_group | unsigned int | 1 | the number of iterations in a DIIS group | DIIS extrapolation is only used for the first \c diis_num_extrap_group of these iterations |
   * | diis_num_extrap_group | unsigned int | 1 | the number of DIIS extrapolations to do at the beginning of an iteration group |
   *
   * example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "wfn": {
   *    "type": "zRHF",
   *    "atoms": "$:h2o",
   *    "wfn_world": "$:wfn_world",
   *    "max_iter": 100,
   *    "soad_guess": true,
   *    "print_detail": true,
   *    "max_condition_num": 1e8,
   *    "print_max_item": 100,
   *    "fock_mixing": 0.0,
   *    "diis": "gamma_point",
   *    "k_points": [1, 1, 11]
   *  }
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  // clang-format on
  zRHF(const KeyVal& kv);

  ~zRHF() {}

  void obsolete() override;

  /// return crystal orbital coefficients
  MatrixzVec co_coeff() override { return C_; }

  /// return crystal orbital energies
  VectordVec co_energy() override { return eps_; }

  /// return # of k points in each direction
  Vector3i nk() override { return nk_; }

  /// return the cardinal number of k points
  int64_t k_size() override { return k_size_; }

 private:
  /*!
   * \brief This performs SCF procedure for zRHF
   *
   * @param thresh the target SCF convergence threshold
   * @todo document \c thresh parameter
   * @throws MaxIterExceeded if exceeded the limit on the number of iteration
   */
  void solve(double thresh);

  /*!
   * \brief This diagonalizes Fock matrix in reciprocal space and
   * computes density: D_ = Int_k( Exp(I k.R) C(occ).C(occ)t )
   */
  std::pair<array_type, array_type_z> compute_density();

  /*!
   * \brief This transforms an integral matrix from real to reciprocal space
   * via M(μ, ν_k) = \sum_R exp(I k.R) M(μ, ν_R).
   * \param matrix the real-space integral matrix M(μ, ν_R)
   * \param real_lattice_range maximum unit cell index (n1, n2, n3) in lattice
   * sum of \c R. n1, n2 and n3 are non-negative integers.
   * \param recip_lattice_range number of k points (k1, k2, k3) in reciprocal
   * space. k1, k2, and k3 are positive integers
   * \return the reciprocal-space integral matrix M(μ, ν_k)
   */
  array_type_z transform_real2recip(const array_type& matrix,
                                    const Vector3i& real_lattice_range,
                                    const Vector3i& recip_lattice_range);

  /*!
   * \brief This transforms an integral matrix from real to reciprocal space
   * via M(μ, ν_k) = \sum_R exp(I k.R) M(μ, ν_R).
   * \param matrix the real-space integral matrix M(μ, ν_R)
   * \return the reciprocal-space integral matrix M(μ, ν_k)
   */
  array_type_z transform_real2recip(const array_type& matrix);

  /*!
   * \brief This changes phase factor of a complex value
   * \param arg_value original complex value
   * \param factor \phi in e^(i \phi)
   */
  MatrixZ reverse_phase_factor(MatrixZ& mat0);

 protected:
  array_type S_;
  array_type D_;
  bool print_detail_;
  int64_t print_max_item_;
  std::unique_ptr<scf::PeriodicFockBuilder<Tile, Policy>> f_builder_;

 private:
  array_type T_;
  array_type V_;
  array_type H_;
  array_type J_;
  array_type K_;
  array_type F_;
  array_type_z Sk_;
  array_type_z Fk_;
  array_type_z Dk_;

  MatrixzVec C_;
  VectordVec eps_;
  MatrixzVec X_;

  double energy_;
  int64_t docc_;

  const KeyVal kv_;
  int64_t maxiter_;
  double max_condition_num_;
  double fmix_;

  std::string diis_;
  unsigned int diis_start_;
  unsigned int diis_num_vecs_;
  double diis_damping_;
  double diis_mixing_;
  unsigned int diis_num_iters_group_;
  unsigned int diis_num_extrap_group_;

  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  Vector3i nk_;
  Vector3d dcell_;
  int64_t R_size_;
  int64_t RJ_size_;
  int64_t RD_size_;
  int64_t k_size_;

  double init_duration_ = 0.0;
  double j_duration_ = 0.0;
  double k_duration_ = 0.0;
  double trans_duration_ = 0.0;
  double d_duration_ = 0.0;
  double scf_duration_ = 0.0;

  /*!
   * \brief This initializes zRHF by assigning values to private members
   * and computing initial guess for the density
   *
   * \param kv KeyVal object
   */
  virtual void init(const KeyVal& kv);

  bool can_evaluate(Energy* energy) override;
  void evaluate(Energy* result) override;

  /// returns Hartree-Fock energy
  virtual double compute_energy();
  /// initializes periodic four-center Fock builder
  virtual void init_fock_builder();

  /// builds Fock
  void build_F();
};

/*!
 * \brief DFzRHF class uses density fitting for Coulomb and 4-center-K for
 * exchange.
 *
 * Refs: Burow, A. M.; Sierka, M.; Mohamed, F. JCP. 131, 214101 (2009)
 */
template <typename Tile, typename Policy>
class DFzRHF : public zRHF<Tile, Policy> {
 public:
  using factory_type = typename zRHF<Tile, Policy>::factory_type;

  /**
   * \brief KeyVal constructor for DFzRHF
   *
   * \param kv The KeyVal object takes same keywords in zRHF.
   *
   * example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "wfn": {
   *    "type": "DF-zRHF",
   *    "atoms": "$:h2o",
   *    "wfn_world": "$:wfn_world",
   *    "max_iter": 100,
   *    "soad_guess": true,
   *    "print_detail": true,
   *    "max_condition_num": 1e8,
   *    "print_max_item": 100,
   *    "fock_mixing": 0.0,
   *    "diis": "gamma_point",
   *    "k_points": [1, 1, 11]
   *  }
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  DFzRHF(const KeyVal& kv);

  ~DFzRHF() {}

 private:
  /// initializes necessary arrays for DFzRHF Fock builder
  void init_fock_builder() override;
};

/*!
 * \brief four-center zRHF class uses shell-level&screening 4-center Fock
 * builder
 */
template <typename Tile, typename Policy>
class FourCenterzRHF : public zRHF<Tile, Policy> {
 public:

  /**
   * \brief KeyVal constructor for FourCenterzRHF
   *
   * \param kv The KeyVal object takes same keywords in zRHF.
   *
   * example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "wfn": {
   *    "type": "FourCenter-zRHF",
   *    "atoms": "$:h2o",
   *    "wfn_world": "$:wfn_world",
   *    "max_iter": 100,
   *    "soad_guess": true,
   *    "print_detail": true,
   *    "max_condition_num": 1e8,
   *    "print_max_item": 100,
   *    "fock_mixing": 0.0,
   *    "diis": "gamma_point",
   *    "k_points": [1, 1, 11]
   *  }
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  FourCenterzRHF(const KeyVal& kv);

  ~FourCenterzRHF() {}

 private:
  void init_fock_builder() override;
};

/*!
 * \brief RIJCADFKzRHF uses RI-J for coulomb and CADF-K for exchange
 */
template <typename Tile, typename Policy>
class RIJCADFKzRHF : public zRHF<Tile, Policy> {
 public:
  using factory_type = typename zRHF<Tile, Policy>::factory_type;

  /**
   * \brief KeyVal constructor for RIJCADFKzRHF
   *
   * \param kv The KeyVal object takes same keywords in zRHF, and the following keywords:
   *  | Keyword | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |\c force_shape_threshold | double | 0.0 | This gives the threshold used to construct the shape of F(Υ, μ, ν) using the shape of Q(Y, ρ, ν). See periodic_cadf_k_builder.h for more details.|
   *
   * example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "wfn": {
   *    "type": "RIJ-CADFK-zRHF",
   *    "atoms": "$:h2o",
   *    "wfn_world": "$:wfn_world",
   *    "max_iter": 100,
   *    "soad_guess": true,
   *    "print_detail": true,
   *    "max_condition_num": 1e8,
   *    "print_max_item": 100,
   *    "force_shape_threshold": 1e-10,
   *    "fock_mixing": 0.0,
   *    "diis": "gamma_point",
   *    "k_points": [1, 1, 11]
   *  }
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  RIJCADFKzRHF(const KeyVal& kv);

  ~RIJCADFKzRHF() {}

 private:
  void init_fock_builder() override;
  double force_shape_threshold_;
};

/*!
 * \brief FourCenterJCADFKzRHF uses four-center-J for coulomb
 * and CADF-K for exchange
 */
template <typename Tile, typename Policy>
class FourCenterJCADFKzRHF : public zRHF<Tile, Policy> {
 public:
  using factory_type = typename zRHF<Tile, Policy>::factory_type;

  /**
   * \brief KeyVal constructor for FourCenterJCADFKzRHF
   *
   * \param kv The KeyVal object takes same keywords in zRHF, and the following keywords:
   *  | Keyword | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |\c force_shape_threshold | double | 0.0 | This gives the threshold used to construct the shape of F(Υ, μ, ν) using the shape of Q(Y, ρ, ν). See periodic_cadf_k_builder.h for more details.|
   *
   * example input:
   *
   * ~~~~~~~~~~~~~~~~~~~~~{.json}
   *  "wfn": {
   *    "type": "FourCenterJ-CADFK-zRHF",
   *    "atoms": "$:h2o",
   *    "wfn_world": "$:wfn_world",
   *    "max_iter": 100,
   *    "soad_guess": true,
   *    "print_detail": true,
   *    "max_condition_num": 1e8,
   *    "print_max_item": 100,
   *    "force_shape_threshold": 1e-10,
   *    "fock_mixing": 0.0,
   *    "diis": "gamma_point",
   *    "k_points": [1, 1, 11]
   *  }
   * ~~~~~~~~~~~~~~~~~~~~~
   */
  FourCenterJCADFKzRHF(const KeyVal& kv);

  ~FourCenterJCADFKzRHF() {}

 private:
  void init_fock_builder() override;
  double force_shape_threshold_;
};

#if TA_DEFAULT_POLICY == 0

#elif TA_DEFAULT_POLICY == 1
extern template class zRHF<TA::TensorD, TA::SparsePolicy>;
extern template class DFzRHF<TA::TensorD, TA::SparsePolicy>;
extern template class FourCenterzRHF<TA::TensorD, TA::SparsePolicy>;
extern template class RIJCADFKzRHF<TA::TensorD, TA::SparsePolicy>;
extern template class FourCenterJCADFKzRHF<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace  lcao
}  // namespace  mpqc

#include "zrhf_impl.h"
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_
