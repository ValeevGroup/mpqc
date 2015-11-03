#pragma once

#ifndef MPQC_INTEGRALS_SCHWARZSCREEN_H
#define MPQC_INTEGRALS_SCHWARZSCREEN_H

#include "../task_integrals_common.h"
#include "../../common/typedefs.h"
#include "screen_base.h"
#include "../../include/tiledarray.h"

#include <stdexcept>
#include "../../include/libint.h"


namespace mpqc {
namespace integrals {

/*! \brief Class which holds shell set information for screening.  */
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
    Qmatrix(MatrixD Q);

    double const &Q(int64_t row, int64_t col) const { return Q_(row, col); }

    double const &max_in_row(int64_t row) const {
        return max_elem_in_row_(row);
    }

    double max_in_Q() const { return max_elem_in_Q_; }
};

/*! \brief Class for Schwarz based screening.
 *
 * We will assume that these screeners are replicated for the integrals they
 * need and so will not bother with serialization of them.
 *
 * We need to hold two Q matrices because sometimes we may want to screen from
 * different bases in the bra and ket, one example being Density Fitting.
 */
class SchwarzScreen : public Screener {
  private:
    std::shared_ptr<Qmatrix> Qab_;
    std::shared_ptr<Qmatrix> Qcd_;
    double thresh_;
    static constexpr int64_t order3_col_index = 0;

  public:
    const double &Qab(int64_t ordA) const {
        assert(Qab_ != nullptr);
        return Qab_->Q(ordA, order3_col_index);
    }

    const double &Qab(int64_t ordA, int64_t ordB) const {
        assert(Qab_ != nullptr);
        return Qab_->Q(ordA, ordB);
    }

    const double &Qcd(int64_t ordC, int64_t ordD) const {
        assert(Qcd_ != nullptr);
        return Qcd_->Q(ordC, ordD);
    }

    double max_cd() const {
        assert(Qcd_ != nullptr);
        return Qcd_->max_in_Q();
    }

    double max_ab() const {
        assert(Qab_ != nullptr);
        return Qab_->max_in_Q();
    }

    double max_in_row_ab(int64_t ord) const {
        assert(Qab_ != nullptr);
        return Qab_->max_in_row(ord);
    }

    double max_in_row_cd(int64_t ord) const {
        assert(Qcd_ != nullptr);
        return Qcd_->max_in_row(ord);
    }

    double four_center_skip(int64_t ordA) const {
        if (Qcd_ != nullptr) {
            return max_in_row_ab(ordA) * max_cd();
        } else {
            return max_in_row_ab(ordA) * max_ab();
        }
    }

    double four_center_skip(int64_t ordA, int64_t ordB) const {
        if (Qcd_ != nullptr) {
            return Qab(ordA, ordB) * max_cd();
        } else {
            return Qab(ordA, ordB) * max_ab();
        }
    }

    double four_center_skip(int64_t ordA, int64_t ordB, int64_t ordC) const {
        if (Qcd_ != nullptr) {
            return Qab(ordA, ordB) * max_in_row_cd(ordC);
        } else {
            return Qab(ordA, ordB) * max_in_row_ab(ordC);
        }
    }

    double four_center_skip(int64_t ordA, int64_t ordB, int64_t ordC,
                            int64_t ordD) const {
        if (Qcd_ != nullptr) {
            return Qab(ordA, ordB) * Qcd(ordC, ordD);
        } else {
            return Qab(ordA, ordB) * Qab(ordC, ordD);
        }
    }

    SchwarzScreen() = default;
    SchwarzScreen(SchwarzScreen const &) = default;
    SchwarzScreen(SchwarzScreen &&) = default;
    SchwarzScreen &operator=(SchwarzScreen const &) = default;
    SchwarzScreen &operator=(SchwarzScreen &&) = default;

    virtual ~SchwarzScreen() = default;

    /*! \brief Constructor which requires a Q matrix.
     *
     * The threshold for screening defaults to 1e-10, but is settable,
     * The second array should be used in the event that the Ket and Bra
     * bases are different.
     *
     * For Density Fitting integrals of type (A|cd) Qab_ represents A and
     * Qcd_ represents c and d, The reverse ordering is not supported.
     */
    SchwarzScreen(std::shared_ptr<Qmatrix> Qab,
                  std::shared_ptr<Qmatrix> Qcd = nullptr, double thresh = 1e-10)
            : Screener(),
              Qab_(std::move(Qab)),
              Qcd_(std::move(Qcd)),
              thresh_(thresh) {}

    double skip_threshold() const { return thresh_; }

    // Not overriding the two index parts

    // Assume for all three center that ordering is (A|cd)
    /// Three loop Outer Screen.
    bool skip(int64_t ordA, Shell const &, ShellVec const &,
              ShellVec const &) override {
        const auto est = Qab(ordA) * max_cd();
        return (est < thresh_) ? true : false;
    }

    /// Three loop Middle Screen
    bool skip(int64_t ordA, int64_t ordC, Shell const &, Shell const &,
              ShellVec const &) override {
        const auto est = Qab(ordA) * max_in_row_cd(ordC);
        return (est < thresh_) ? true : false;
    }

    /// Three loop Inner Screen.
    bool skip(int64_t ordA, int64_t ordC, int64_t ordD, Shell const &,
              Shell const &, Shell const &) override {
        const auto est = Qab(ordA) * Qcd(ordC, ordD);
        return (est < thresh_) ? true : false;
    }

    /// Four loop Outer Screen.
    bool skip(int64_t ordA, Shell const &, ShellVec const &, ShellVec const &,
              ShellVec const &) override {
        const auto est = four_center_skip(ordA);
        return (est < thresh_) ? true : false;
    }

    /// Four loop Second Screen.
    bool skip(int64_t ordA, int64_t ordB, Shell const &, Shell const &,
              ShellVec const &, ShellVec const &) override {
        const auto est = four_center_skip(ordA, ordB);
        return (est < thresh_) ? true : false;
    }

    /// Four loop Third Screen.
    bool skip(int64_t ordA, int64_t ordB, int64_t ordC, Shell const &,
              Shell const &, Shell const &, ShellVec const &) override {
        const auto est = four_center_skip(ordA, ordB, ordC);
        return (est < thresh_) ? true : false;
    }

    /// Four loop Inner Screen.
    bool
    skip(int64_t ordA, int64_t ordB, int64_t ordC, int64_t ordD, Shell const &,
         Shell const &, Shell const &, Shell const &) override {
        const auto est = four_center_skip(ordA, ordB, ordC, ordD);
        return (est < thresh_) ? true : false;
    }
};

namespace detail {

/*! \brief Function for making Q for Density Fitting auxiliary basis.
 *
 * Returns a Vector and thus should not be used for basis sets with two indices.
 *
 */
template <typename E>
inline std::shared_ptr<Qmatrix>
compute_Q(madness::World &world, ShrPool<E> const &eng,
          basis::Basis const &bs0) {

    const auto shells = bs0.flattened_shells();
    const auto nshells = shells.size();

    VectorD Q(nshells);

    // pass by pointers since tasks copy params
    auto task_func = [=](Shell const *sh, double *Q_val) {
        const auto *buf
              = eng->local().compute(*sh, unit_shell, *sh, unit_shell);

        const auto nsh = sh->size();
        const auto bmap = Eig::Map<const MatrixD>(buf, nsh, nsh);

        *Q_val = std::sqrt(bmap.lpNorm<2>());
    };

    for (auto i = 0ul; i < nshells; ++i) {
        auto *Q_val_ptr = &Q(i);
        world.taskq.add(task_func, &shells[i], Q_val_ptr);
    }
    world.gop.fence();

    return std::make_shared<Qmatrix>(std::move(Q));
}

// Just a type to stand in as a parameter
struct form_matrix_dummy {};

/*! \brief Function for making Q for a two index basis
 *
 *  Use dummy parameter to select this overload.
 */
template <typename E>
inline std::shared_ptr<Qmatrix>
compute_Q(madness::World &world, ShrPool<E> const &eng, basis::Basis const &bs,
          form_matrix_dummy const &) {

    auto const shells0 = bs.flattened_shells();
    const auto nshells0 = shells0.size();

    MatrixD Q(nshells0, nshells0);

    // pass by pointers since tasks copy params
    auto task_func = [=](Shell const *sh0, Shell const *sh1, double *Q_val) {
        const auto *buf = eng->local().compute(*sh0, *sh1, *sh0, *sh1);

        const auto n2 = sh0->size() * sh1->size();
        const auto bmap = Eig::Map<const MatrixD>(buf, n2, n2);

        *Q_val = std::sqrt(bmap.lpNorm<2>());
    };

    for (auto i = 0ul; i < nshells0; ++i) {
        for (auto j = 0ul; j < nshells0; ++j) {
            auto *Q_val_ptr = &Q(i, j);
            world.taskq.add(task_func, &shells0[i], &shells0[j], Q_val_ptr);
        }
    }
    world.gop.fence();

    return std::make_shared<Qmatrix>(std::move(Q));
}

/*! \brief Function for making Q for a two index basis */
template <typename E>
inline std::shared_ptr<Qmatrix>
compute_Q(madness::World &world, ShrPool<E> const &eng, basis::Basis const &bs0,
          basis::Basis const &bs1) {

    auto const shells0 = bs0.flattened_shells();
    auto const shells1 = bs1.flattened_shells();
    const auto nshells0 = shells0.size();
    const auto nshells1 = shells1.size();

    MatrixD Q(nshells0, nshells1);


    // pass by pointers since tasks copy params
    auto task_func = [=](Shell const *sh0, Shell const *sh1, double *Q_val) {
        const auto *buf = eng->local().compute(*sh0, *sh1, *sh0, *sh1);

        const auto n2 = sh0->size() * sh1->size();
        const auto bmap = Eig::Map<const MatrixD>(buf, n2, n2);

        *Q_val = std::sqrt(bmap.lpNorm<2>());
    };

    for (auto i = 0ul; i < nshells0; ++i) {
        for (auto j = 0ul; j < nshells1; ++j) {
            auto *Q_val_ptr = &Q(i, j);
            world.taskq.add(task_func, &shells0[i], &shells1[j], Q_val_ptr);
        }
    }
    world.gop.fence();

    return std::make_shared<Qmatrix>(std::move(Q));
}


} // namespace detail;


/*! \brief struct which builds SchwarzScreen screeners */
class init_schwarz_screen {
  private:
    double threshold = 1e-10;

    // Add more overloads as desired, can have alternate ways of passing bases.
    template <typename E>
    SchwarzScreen
    compute_df(madness::World &world, ShrPool<E> &eng, basis::Basis const &dfbs,
               basis::Basis const &obs) const {
        auto Q_a = detail::compute_Q(world, eng, dfbs);
        auto Q_cd
              = detail::compute_Q(world, eng, obs, detail::form_matrix_dummy{});

        return SchwarzScreen(std::move(Q_a), std::move(Q_cd), threshold);
    }

    // This guy is here for interface reasons only. 
    template <typename E>
    SchwarzScreen
    compute_df(madness::World &, ShrPool<E> &, basis::Basis const &) const {
        assert(false);
        return SchwarzScreen();
    }

    template <typename E>
    SchwarzScreen
    compute_4c(madness::World &world, ShrPool<E> &eng, basis::Basis const bs) const {
        auto Q_ab
              = detail::compute_Q(world, eng, bs, detail::form_matrix_dummy{});
        return SchwarzScreen(std::move(Q_ab), nullptr, threshold);
    }

  public:
    init_schwarz_screen() = default;
    init_schwarz_screen(double thresh) : threshold(thresh) {}

    enum class ScreenType { DensityFitting, FourCenter };

    template <typename E, typename... Bases>
    SchwarzScreen operator()(madness::World &world, ShrPool<E> &eng,
                             ScreenType const &type, Bases &&... bases) const {
        if (type == ScreenType::DensityFitting) {
            // need at least 2 basis sets for DF
            assert(sizeof...(Bases) >= 2);
            return compute_df(world, eng, std::forward<Bases>(bases)...);
        }
        else {
            return compute_4c(world, eng, std::forward<Bases>(bases)...);
        }
    }
};


} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_SCHWARZSCREEN_H
