#pragma once

#ifndef MPQC_INTEGRALS_SCHWARZSCREEN_H
#define MPQC_INTEGRALS_SCHWARZSCREEN_H

#include "task_integrals_common.h"
#include "../common/typedefs.h"
#include "integral_screeners.h"
#include "../include/tiledarray.h"

#include "../include/libint.h"
#include <stdexcept>

namespace mpqc {
namespace integrals {


class Qmatrix {
  private:
    MatrixD Q_;
    VectorD max_elem_in_row_;
    double max_elem_in_Q_;

  public:
    Qmatrix() = default;
    Qmatrix(Qmatrix const &) = default;
    Qmatrix(Qmatrix &&) = default;
    Qmatrix &operator=(Qmatrix const &) = default;
    Qmatrix &operator=(Qmatrix &&) = default;
    Qmatrix(MatrixD Q)
            : Q_(std::move(Q)),
              max_elem_in_row_(VectorD(Q_.rows())),
              max_elem_in_Q_(0.0)

    {
        const auto nrows = Q_.rows();
        for (auto i = 0; i < nrows; ++i) {
            const auto max_elem = Q_.row(i).cwiseAbs().maxCoeff();
            max_elem_in_Q_ = std::max(max_elem, max_elem_in_Q_);
            max_elem_in_row_(i) = max_elem;
        }
    }

    double const &Q(int64_t row, int64_t col) const { return Q_(row, col); }

    double const &max_in_row(int64_t row) const {
        return max_elem_in_row_(row);
    }

    double max_in_Q() const { return max_elem_in_Q_; }
};

class SchwarzScreen : public Screener {
  private:
    std::shared_ptr<Qmatrix> Qab_;
    std::shared_ptr<Qmatrix> Qcd_;
    double thresh_;
    static constexpr int64_t order3_col_index = 0;

  public:
    SchwarzScreen() = default;
    SchwarzScreen(SchwarzScreen const &) = default;
    SchwarzScreen(SchwarzScreen &&) = default;
    SchwarzScreen &operator=(SchwarzScreen const &) = default;
    SchwarzScreen &operator=(SchwarzScreen &&) = default;

    virtual ~SchwarzScreen() = default;

    SchwarzScreen(std::shared_ptr<Qmatrix> Qab, std::shared_ptr<Qmatrix> Qcd,
                  double thresh)
            : Screener(),
              Qab_(std::move(Qab)),
              Qcd_(std::move(Qcd)),
              thresh_(thresh) {}

    // Not overriding the two index parts

    // Assume for all three center that ordering is (A|cd)
    /// Three loop Outer Screen.
    bool skip(int64_t ordA, Shell const &, ShellVec const &,
              ShellVec const &) override {
        const auto est = Qab_->Q(ordA, order3_col_index) * Qcd_->max_in_Q();
        return (est < thresh_) ? true : false;
    }

    /// Three loop Middle Screen
    bool skip(int64_t ordA, int64_t ordC, Shell const &, Shell const &,
              ShellVec const &) override {
        const auto est = Qab_->Q(ordA, order3_col_index)
                         * Qcd_->max_in_row(ordC);
        return (est < thresh_) ? true : false;
    }

    /// Three loop Inner Screen.
    bool skip(int64_t ordA, int64_t ordC, int64_t ordD, Shell const &,
              Shell const &, Shell const &) override {
        const auto est = Qab_->Q(ordA, order3_col_index) * Qcd_->Q(ordC, ordD);
        return (est < thresh_) ? true : false;
    }

    /// Four loop Outer Screen.
    bool skip(int64_t ordA, Shell const &, ShellVec const &, ShellVec const &,
              ShellVec const &) override {
        const auto est = Qab_->max_in_row(ordA) * Qcd_->max_in_Q();
        return (est < thresh_) ? true : false;
    }

    /// Four loop Second Screen.
    bool skip(int64_t ordA, int64_t ordB, Shell const &, Shell const &,
              ShellVec const &, ShellVec const &) override {
        const auto est = Qab_->Q(ordA, ordB) * Qcd_->max_in_Q();
        return (est < thresh_) ? true : false;
    }

    /// Four loop Third Screen.
    bool skip(int64_t ordA, int64_t ordB, int64_t ordC, Shell const &,
              Shell const &, Shell const &, ShellVec const &) override {

        const auto est = Qab_->Q(ordA, ordB) * Qcd_->max_in_row(ordC);
        return (est < thresh_) ? true : false;
    }

    /// Four loop Inner Screen.
    bool
    skip(int64_t ordA, int64_t ordB, int64_t ordC, int64_t ordD, Shell const &,
         Shell const &, Shell const &, Shell const &) override {
        const auto est = Qab_->Q(ordA, ordB) * Qcd_->Q(ordC, ordD);
        return (est < thresh_) ? true : false;
    }
};

namespace detail {

template <typename E>
inline std::shared_ptr<Qmatrix>
compute_Q(ShrPool<E> const &eng, int64_t index, basis::Basis const &bs0) {
    auto const &shellvec = bs0.cluster_shells()[index];
    const auto nshells = shellvec.size();

    VectorD Q(nshells);
    for (auto i = 0ul; i < nshells; ++i) {
        auto const &sh = shellvec[i];
        const auto nsh = sh.size();

        const auto *buf = eng->local().compute(sh, unit_shell, sh, unit_shell);

        const auto bmap = Eig::Map<const MatrixD>(buf, nsh, nsh);
        Q(i) = std::sqrt(bmap.lpNorm<2>());
    }

    return std::make_shared<Qmatrix>(std::move(Q));
}

template <typename E>
inline std::shared_ptr<Qmatrix>
compute_Q(ShrPool<E> const &eng, int64_t index0, int64_t index1,
          basis::Basis const &bs0, basis::Basis const &bs1) {
    auto const &shellvec0 = bs0.cluster_shells()[index0];
    auto const &shellvec1 = bs1.cluster_shells()[index1];
    const auto nshells0 = shellvec0.size();
    const auto nshells1 = shellvec1.size();

    MatrixD Q(nshells0, nshells1);

    for (auto i = 0ul; i < nshells0; ++i) {
        auto const &sh0 = shellvec0[i];
        const auto nsh0 = sh0.size();

        for (auto j = 0ul; j < nshells1; ++j) {
            auto const &sh1 = shellvec1[j];
            const auto nsh1 = sh1.size();
            const auto n2 = nsh0 * nsh1;

            const auto *buf = eng->local().compute(sh0, sh1, sh0, sh1);

            const auto bmap = Eig::Map<const MatrixD>(buf, n2, n2);
            Q(i, j) = std::sqrt(bmap.lpNorm<2>());
        }
    }

    return std::make_shared<Qmatrix>(std::move(Q));
}

} // namespace detail;


struct init_schwarz_screen {
    template <typename E>
    SchwarzScreen operator()(detail::IdxVec const &,
                             detail::ShrBases<2> const &, ShrPool<E> const &) {
        throw std::logic_error("Cannot use SchwarzScreen class on Tensors with "
                               "only two dimensions.");
        return SchwarzScreen();
    }

    template <typename E>
    SchwarzScreen
    operator()(detail::IdxVec const &idx,
               detail::ShrBases<3> const &bases, ShrPool<E> const &engs) {
        auto Q_a = detail::compute_Q(engs, idx[0], bases->operator[](0));
        auto Q_cd
              = detail::compute_Q(engs, idx[1], idx[2], bases->operator[](1),
                                  bases->operator[](2));

        auto thresh = SpShapeF::threshold();

        return SchwarzScreen(std::move(Q_a), std::move(Q_cd), thresh);
    }

    template <typename E>
    SchwarzScreen
    operator()(detail::IdxVec const &idx,
               detail::ShrBases<4> const &bases, ShrPool<E> const &engs) {
        auto Q_ab
              = detail::compute_Q(engs, idx[0], idx[1], bases->operator[](0),
                                  bases->operator[](1));

        auto Q_cd
              = detail::compute_Q(engs, idx[2], idx[3], bases->operator[](2),
                                  bases->operator[](3));

        auto thresh = SpShapeF::threshold();

        return SchwarzScreen(std::move(Q_ab), std::move(Q_cd), thresh);
    }
};


} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_SCHWARZSCREEN_H
