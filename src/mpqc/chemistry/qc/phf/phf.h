#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PHF_PHF_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PHF_PHF_H_

#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

#include <memory>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/wfn/ao_wfn.h"
#include "mpqc/util/external/c++/memory"

namespace mpqc {
namespace phf {

/**
 * Periodic Hartree-Fock class
 */

class PHF : public qc::Wavefunction {
public:
    using TArray = TA::DistArray<TA::TensorZ, TA::SparsePolicy>;
    using PeriodicAOIntegral = integrals::PeriodicAOFactory<TA::TensorZ, TA::SparsePolicy>;
    using MatrixcVec = std::vector<Matrixc>;

    PHF() = default;

    /**
     * KeyVal constructor for PHF
     *
     * keywords: takes all keywords from AOWavefunction
     *
     * | KeyWord | Type | Default| Description |
     * |---------|------|--------|-------------|
     * | converge | double | 1.0e-07 | converge limit |
     * | max_iter | int | 30 | maximum number of iteration |
     *
     */
    PHF(const KeyVal& kv);

    ~PHF() = default;

    void compute(qc::PropertyBase *pb) override;
    void obsolete() override;

    double value() override;

    bool solve();

private:
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
    TArray C_;

    MatrixcVec X_;

    double repulsion_;
    int64_t docc_;

    const KeyVal kv_;
    double converge_;
    int64_t maxiter_;
    bool print_detail_;

    double init_duration_ = 0.0;
    double j_duration_ = 0.0;
    double k_duration_ = 0.0;
    double trans_duration_ = 0.0;
    double d_duration_ = 0.0;
    double scf_duration_ = 0.0;

    void init(const KeyVal& kv);
};

} // end of namespace phf
} // end of namespace mpqc
#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_PHF_PHF_H_
