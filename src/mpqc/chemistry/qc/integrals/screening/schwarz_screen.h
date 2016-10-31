#pragma once

#ifndef MPQC_INTEGRALS_SCHWARZSCREEN_H
#define MPQC_INTEGRALS_SCHWARZSCREEN_H

#include <tiledarray.h>

#include "../../../../../../common/typedefs.h"
#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>
#include <mpqc/chemistry/qc/integrals/screening/screen_base.h>

#include <stdexcept>
#include <unordered_map>


namespace mpqc {
namespace integrals {

/*! \brief Class which holds shell set information for screening.
 *
 * Must keep a function to shell map which takes the index of a function and
 * returns the index of its shell.  The map only maps the first function in
 * each shell so should fail if say the global index of the second function in
 * the shell is passed in.
 */
class Qmatrix {
  private:
    MatrixD Q_;
    std::unordered_map<int64_t, int64_t> func_to_shell_map_;
    VectorD max_elem_in_row_;
    double max_elem_in_Q_;

    // Convert a function index into a shell index
    int64_t f2s(int64_t a) const {
        auto it = func_to_shell_map_.find(a);

        // If this hits the issue was likely an index that was not the
        // first function in the shell.
        assert(it != func_to_shell_map_.end());
        return it->second;
    }

  public:
    Qmatrix() = default;
    Qmatrix(Qmatrix const &) = default;
    Qmatrix(Qmatrix &&) = default;
    Qmatrix &operator=(Qmatrix const &) = default;
    Qmatrix &operator=(Qmatrix &&) = default;
    Qmatrix(MatrixD Q, std::unordered_map<int64_t, int64_t>);

    // Get whole matrix
    MatrixD const &Q() const { return Q_; }
    // If vector get index
    double const &Q(int64_t idx) const { return Q_(f2s(idx)); }
    // If Matrix get elem
    double const &Q(int64_t row, int64_t col) const {
        return Q_(f2s(row), f2s(col));
    }

    double const &max_in_row(int64_t row) const {
        return max_elem_in_row_(f2s(row));
    }

    double max_in_Q() const { return max_elem_in_Q_; }

    int64_t func_to_shell(int64_t func_idx) const { return f2s(func_idx); }
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

  public:
    bool Qab_is_vector() const { return Qab_->Q().cols() == 1; }

    const double &Qab(int64_t a) const {
        assert(Qab_ != nullptr && Qab_->Q().cols() == 1);
        return Qab_->Q(a);
    }

    const double &Qcd(int64_t a) const {
        assert(Qcd_ != nullptr && Qcd_->Q().cols() == 1);
        return Qab_->Q(a);
    }

    const double &Qab(int64_t a, int64_t b) const {
        assert(Qab_ != nullptr);
        return Qab_->Q(a, b);
    }

    const double &Qcd(int64_t c, int64_t d) const {
        assert(Qcd_ != nullptr);
        return Qcd_->Q(c, d);
    }

    double max_ab() const {
        assert(Qab_ != nullptr);
        return Qab_->max_in_Q();
    }

    double max_cd() const {
        assert(Qcd_ != nullptr);
        return Qcd_->max_in_Q();
    }


    double max_in_row_ab(int64_t row_a) const {
        assert(Qab_ != nullptr);
        return Qab_->max_in_row(row_a);
    }

    double max_in_row_cd(int64_t row_a) const {
        assert(Qcd_ != nullptr);
        return Qcd_->max_in_row(row_a);
    }

    // Screen three center based on first index
    double three_center_est(int64_t a) const { return Qab(a) * max_cd(); }

    double three_center_est(int64_t a, int64_t b) const {
        return Qab(a) * max_in_row_cd(b);
    }

    double three_center_est(int64_t a, int64_t b, int64_t c) const {
        return Qab(a) * Qcd(b, c);
    }

    // Here for interface reasons not actually usable
    double three_center_est(int64_t, int64_t, int64_t, int64_t) const {
        assert(false);
        return 0.0;
    }

    double four_center_est(int64_t a) const {
        return max_in_row_ab(a) * max_ab();
    }

    double four_center_est(int64_t a, int64_t b) const {
        return Qab(a, b) * max_ab();
    }

    double four_center_est(int64_t a, int64_t b, int64_t c) const {
        return Qab(a, b) * max_in_row_ab(c);
    }

    double four_center_est(int64_t a, int64_t b, int64_t c, int64_t d) const {
        return Qab(a, b) * Qab(c, d);
    }

    double four_center_est_2Q(int64_t a) const {
        return max_in_row_ab(a) * max_cd();
    }

    double four_center_est_2Q(int64_t a, int64_t b) const {
        return Qab(a, b) * max_cd();
    }

    double four_center_est_2Q(int64_t a, int64_t b, int64_t c) const {
        return Qab(a, b) * max_in_row_cd(c);
    }

    double
    four_center_est_2Q(int64_t a, int64_t b, int64_t c, int64_t d) const {
        return Qab(a, b) * Qcd(c, d);
    }

    /// Returns the shell index given the starting funciton index of that shell.
    int64_t Qab_f2s(int64_t func_idx) const {
        return Qab_->func_to_shell(func_idx);
    }

    /// Returns the shell index given the starting funciton index of that shell.
    int64_t Qcd_f2s(int64_t func_idx) const {
        return Qcd_->func_to_shell(func_idx);
    }


    template <typename... IDX>
    double screen(IDX... idx) const {
        if (Qab_is_vector()) { // Three Center Integrals
            assert(Qcd_ != nullptr);
            assert(sizeof...(IDX) <= 3);
            return three_center_est(idx...);
        } else if (Qcd_ == nullptr) { // Four center with only 1 matrix
            return four_center_est(idx...);
        } else { // Four center with 2 matrices
            return four_center_est_2Q(idx...);
        }
    }

    SchwarzScreen() = default;
    SchwarzScreen(SchwarzScreen const &) = default;
    SchwarzScreen(SchwarzScreen &&) = default;
    SchwarzScreen &operator=(SchwarzScreen const &) = default;
    SchwarzScreen &operator=(SchwarzScreen &&) = default;

    virtual ~SchwarzScreen() noexcept = default;

    /*! \brief Constructor which requires a Q matrix.
     *
     * The threshold for screening defaults to 1e-10, but is settable,
     * The second array should be used in the event that the Ket and Bra
     * bases are different.
     *
     * For Density Fitting integrals of type (A|cd) Qab_ represents A and
     * Qcd_ represents c and d, The reverse ordering is not supported.
     * Qab_ will be a vector.
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
    bool skip(int64_t a) override {
        return (screen(a) < thresh_) ? true : false;
    }

    bool skip(int64_t a, int64_t b) override {
        return (screen(a, b) < thresh_) ? true : false;
    }

    bool skip(int64_t a, int64_t b, int64_t c) override {
        return (screen(a, b, c) < thresh_) ? true : false;
    }

    bool skip(int64_t a, int64_t b, int64_t c, int64_t d) override {
        return (screen(a, b, c, d) < thresh_) ? true : false;
    }

};

namespace detail {

// Make map to shell idx given the first function in the shell.
inline std::unordered_map<int64_t, int64_t>
func_to_shell(std::vector<Shell> const &shells) {
    std::unordered_map<int64_t, int64_t> f2s;
    const auto size = shells.size();

    auto func_idx = 0;
    for (auto i = 0ul; i < size; ++i) {
        f2s.emplace(func_idx, i);
        func_idx += shells[i].size();
    }

    return f2s;
}

// Compute Q vector
template <typename E>
inline MatrixD auxiliary_Q(madness::World &world, ShrPool<E> const &eng,
                           std::vector<Shell> const &shells) {

    // pass by pointers since tasks copy params
    auto task_func = [=](Shell const *sh, Shell const *ush, double *Q_val) {
        const auto& bufs = eng->local().compute(*sh, *ush, *sh, *ush);
        TA_USER_ASSERT(bufs.size() == 1, "unexpected result from Engine::compute");
        const auto nsh = sh->size();
        const auto bmap = Eig::Map<const MatrixD>(bufs[0], nsh, nsh);

        *Q_val = std::sqrt(bmap.lpNorm<2>());
    };

    const auto nshells = shells.size();
    VectorD Q(nshells);

    auto const *ush = &unit_shell;
    for (auto i = 0ul; i < nshells; ++i) {
        auto *Q_val_ptr = &Q(i);
        world.taskq.add(task_func, &shells[i], ush, Q_val_ptr);
    }
    world.gop.fence();

    return Q;
}

// Compute Q matrix
template <typename E>
inline MatrixD four_center_Q(madness::World &world, ShrPool<E> const &eng,
                             std::vector<Shell> const &shells) {

    auto task_func =
          [=](Shell const *sh0, Shell const *sh1, double *Q_val) {
        eng->local().set_precision(0.);
        const auto& bufs = eng->local().compute(*sh0, *sh1, *sh0, *sh1);
        TA_USER_ASSERT(bufs.size() == 1,
                       "unexpected result from Engine::compute");

        const auto n2 = sh0->size() * sh1->size();
        const auto bmap = Eig::Map<const MatrixD>(bufs[0], n2, n2);

        eng->local().set_precision(std::numeric_limits<double>::epsilon());

        *Q_val = std::sqrt(bmap.lpNorm<2>());
    };

    const auto nshells = shells.size();
    MatrixD Q(nshells, nshells);

    for (auto i = 0ul; i < nshells; ++i) {
        for (auto j = 0ul; j < nshells; ++j) {
            auto *Q_val_ptr = &Q(i, j);
            world.taskq.add(task_func, &shells[i], &shells[j], Q_val_ptr);
        }
    }

    world.gop.fence();

    return Q;
}

// Compute Q
template <typename E>
inline std::shared_ptr<Qmatrix>
compute_Q(madness::World &world, ShrPool<E> const &eng, basis::Basis const &bs,
          bool auxiliary_basis = false) {
    const auto shells = bs.flattened_shells();
    MatrixD Q;
    if (auxiliary_basis) {
        Q = auxiliary_Q(world, eng, shells);
    } else {
        Q = four_center_Q(world, eng, shells);
    }

    return std::make_shared<Qmatrix>(std::move(Q), func_to_shell(shells));
}

} // namespace detail;

/*! \brief struct which builds SchwarzScreen screeners */
class init_schwarz_screen {
  private:
    double threshold = 1e-12;

    // Add more overloads as desired, can have alternate ways of passing bases.
    template <typename E>
    SchwarzScreen
    compute_df(madness::World &world, ShrPool<E> &eng, basis::Basis const &dfbs,
               basis::Basis const &obs) const {
        auto Q_a = detail::compute_Q(world, eng, dfbs, true);
        auto Q_cd = detail::compute_Q(world, eng, obs);

        return SchwarzScreen(std::move(Q_a), std::move(Q_cd), threshold);
    }

    template <typename E>
    SchwarzScreen compute_4c(madness::World &world, ShrPool<E> &eng,
                             basis::Basis const &bs) const {
        auto Q_ab = detail::compute_Q(world, eng, bs);
        return SchwarzScreen(std::move(Q_ab), nullptr, threshold);
    }

    // TODO finish this guy mixed basis four center
    template <typename E>
    SchwarzScreen compute_4c(madness::World &, ShrPool<E> &,
                             basis::Basis const &, basis::Basis const &) const {
        assert(false); // Feature not implemented yet.
        return SchwarzScreen();
    }

  public:
    init_schwarz_screen() = default;
    init_schwarz_screen(double thresh) : threshold(thresh) {}

    template <typename E>
    SchwarzScreen operator()(madness::World &world, ShrPool<E> &eng,
                             basis::Basis const &bs) const {
        return compute_4c(world, eng, bs);
    }

    template <typename E>
    SchwarzScreen operator()(madness::World &world, ShrPool<E> &eng,
                             basis::Basis const &bs0, basis::Basis const &bs1,
                             bool mixed_basis_four_center = false) const {
        if (mixed_basis_four_center) {
            assert(false); // This is not yet supported by helper machinary
            return compute_4c(world, eng, bs0, bs1);
        } else {
            return compute_df(world, eng, bs0, bs1);
        }
    }
};


} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_SCHWARZSCREEN_H
