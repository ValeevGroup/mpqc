#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_

#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

#include <memory>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/wfn/ao_wfn.h"
#include "mpqc/util/external/c++/memory"

namespace mpqc {
namespace lcao {

/**
 * complex-valued Restricted Hartree-Fock class
 */

class zRHF : public PeriodicAOWavefunction<TA::TensorZ, TA::SparsePolicy> {
public:
    using Tile = TA::TensorZ;
    using TArray = PeriodicAOWavefunction<TA::TensorZ, TA::SparsePolicy>::ArrayType;
    using PeriodicAOIntegral = PeriodicAOWavefunction::AOIntegral;
    using MatrixcVec = std::vector<Matrixz>;
    using VectorcVec = std::vector<Vectorz>;

    zRHF() = default;

    /**
     * KeyVal constructor for zRHF
     *
     * keywords: takes all keywords from PeriodicAOWavefunction
     *
     * | KeyWord | Type | Default| Description |
     * |---------|------|--------|-------------|
     * | converge | double | 1.0e-07 | converge limit |
     * | max_iter | int | 30 | maximum number of iteration |
     * | soad_guess | bool | true | if use SOAD guess for initial Fock build |
     * | print_detail | bool | false | if print extra computation&time info |
     * | max_condition_num | double | 1.0e8 | maximum condition number for overlap matrix |
     *
     */
    zRHF(const KeyVal& kv);

    ~zRHF() = default;

    void compute(PropertyBase *pb) override;
    void obsolete() override;

    /*!
     * \brief This performs a Hartree-Fock computation
     * \return the Hartree-Fock energy
     */
    double value() override;

private:

    /*!
     * \brief This performs SCF procedure for zRHF
     * \return true if SCF converges, false if not
     */
    bool solve();

    /*!
     * \brief This diagonalizes Fock matrix in reciprocal space and
     * computes density: D_ = Int_k( Exp(I k.R) C(occ).C(occ)t )
     */
    void compute_density();

    TArray T_;
    TArray V_;
    TArray S_;
    TArray Sk_;
    TArray H_;
    TArray Hk_;
    TArray J_;
    TArray K_;
    TArray F_;
    TArray Fk_;
    TArray D_;

    MatrixcVec C_;
    VectorcVec eps_;
    MatrixcVec X_;

    double repulsion_;
    int64_t docc_;

    const KeyVal kv_;
    double converge_;
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
    void init(const KeyVal& kv);
};

}  // namespace  lcao
}  // namespace  mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_ZRHF_H_
