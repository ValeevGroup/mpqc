#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"

#include <memory>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/wfn/ao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/external/c++/memory"

namespace mpqc {
namespace lcao {

using MatrixzVec = std::vector<Matrixz>;
using VectorzVec = std::vector<Vectorz>;
using VectordVec = std::vector<Vectord>;
using Matrix = RowMatrixXd;

/**
 * complex-valued Restricted Hartree-Fock class
 */

class zRHF : public PeriodicAOWavefunction<TA::TensorD, TA::SparsePolicy>,
             public Provides<Energy /*,
          CanonicalOrbitalSpace<TA::DistArray<TA::TensorZ, TA::SparsePolicy>>,
          PopulatedOrbitalSpace<TA::DistArray<TA::TensorZ, TA::SparsePolicy>>*/> {
 public:
  using Tile = TA::TensorD;
  using Policy = TA::SparsePolicy;
  using TArray = PeriodicAOWavefunction<Tile, Policy>::ArrayType;
  using TArrayZ = TA::DistArray<TA::TensorZ, Policy>;
  using factory_type = PeriodicAOWavefunction<Tile, Policy>::AOIntegral;

  zRHF() = default;

  /**
   * KeyVal constructor for zRHF
   *
   * keywords: takes all keywords from PeriodicAOWavefunction
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | converge | double | 1.0e-07 | converge limit |
   * | max_iter | int | 30 | maximum number of iteration |
   * | soad_guess | bool | true | if use SOAD guess for initial Fock build |
   * | print_detail | bool | false | if print extra computation&time info |
   * | max_condition_num | double | 1.0e8 | maximum condition number for overlap
   * matrix |
   *
   */
  zRHF(const KeyVal& kv);

  ~zRHF() = default;

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
  TArray compute_density();

  /*!
   * \brief This transforms an integral matrix from real to reciprocal space
   * \param matrix the real-space integral matrix
   * \return the reciprocal-space integral matrix
   */
  TArrayZ transform_real2recip(TArray& matrix);

  /*!
   * \brief This changes phase factor of a complex value
   * \param arg_value original complex value
   * \param factor \phi in e^(i \phi)
   */
  Matrixz reverse_phase_factor(Matrixz& mat0);

  TArray T_;
  TArray V_;
  TArray S_;
  TArrayZ Sk_;
  TArray H_;
  TArray J_;
  TArray K_;
  TArray F_;
  TArrayZ Fk_;
  TArray D_;

  MatrixzVec C_;
  VectordVec eps_;
  MatrixzVec X_;

  double energy_;
  int64_t docc_;

  const KeyVal kv_;
  int64_t maxiter_;
  bool print_detail_;
  double max_condition_num_;

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
   * \brief This initialize zRHF by assigning values to private members
   * and computing initial guess for the density
   *
   * \param kv KeyVal object
   */
  virtual void init(const KeyVal& kv);

  bool can_evaluate(Energy* energy) override;
  void evaluate(Energy* result) override;

  /// returns Coulomb term J_μν
  virtual TArray J_builder();
  /// returns Exchange term K_μν
  virtual TArray K_builder();
};

/*!
 * \brief DFzRHF class uses density fitting for Coulomb
 *
 * Refs: Burow, A. M.; Sierka, M.; Mohamed, F. JCP. 131, 214101 (2009)
 */
class DFzRHF : public zRHF {
 public:
    using DirectTArray = PeriodicAOWavefunction<Tile, Policy>::DirectTArray;

    DFzRHF(const KeyVal& kv);

    ~DFzRHF() = default;

 private:

    DirectTArray Gamma_;  // (κ λ | G| K) 3-center 2-electron direct integrals
    TArray V_;  // (K |G| Λ) 2-center 2-electron integrals
    TArray P_para_;  // projection matrix that projects X onto auxiliary charge vector
    TArray P_perp_;  // projection matrix that projects X onto the subspace orthogonal to auxiliary charge vector
    TArray C_;  // fitting coefficients
    TArray C_para_;  // the part of C_ that is along with auxiliary charge vector
    TArray C_perp_;  // the part of C_ that is orthogonal to auxiliary charge vector
    TArray M_;  // charge matrix of product density
};

}  // namespace  lcao
}  // namespace  mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_
