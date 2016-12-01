#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PRHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PRHF_H_

#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

#include <memory>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/wfn/ao_wfn.h"
#include "mpqc/util/external/c++/memory"

namespace mpqc {
namespace scf {

/**
 * Periodic Restricted Hartree-Fock class
 */

class PRHF : public qc::Wavefunction {
public:
    using Tile = TA::TensorZ;
    using TArray = TA::DistArray<Tile, TA::SparsePolicy>;
    using PeriodicAOIntegral = integrals::PeriodicAOFactory<TA::TensorZ, TA::SparsePolicy>;
    using MatrixcVec = std::vector<Matrixc>;
    using VectorcVec = std::vector<Vectorc>;

    PRHF() = default;

    /**
     * KeyVal constructor for PRHF
     *
     * keywords: takes all keywords from AOWavefunction
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
    PRHF(const KeyVal& kv);

    ~PRHF() = default;

    void compute(qc::PropertyBase *pb) override;
    void obsolete() override;

    /*!
     * \brief This performs a periodic Hartree-Fock computation
     * \return periodic Hartree-Fock energy
     */
    double value() override;

private:

    /*!
     * \brief This performs SCF procedure for PRHF
     * \return true if SCF converges, false if not
     */
    bool solve();

    /*!
     * \brief This diagonalizes Fock matrix in reciprocal space and
     * computes density: D_ = Int_k( Exp(I k.R) C(occ).C(occ)t )
     */
    void compute_density();

    PeriodicAOIntegral pao_factory_;

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
     * \brief This initialize PRHF by assigning values to private members
     * and computing initial guess for SCF
     *
     * \param kv KeyVal object
     */
    void init(const KeyVal& kv);
};

} // end of namespace scf
} // end of namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PRHF_H_
