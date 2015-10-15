#pragma once

#ifndef MPQC_INTEGRALS_INTEGRALSCREENERS_H
#define MPQC_INTEGRALS_INTEGRALSCREENERS_H

#include "../../common/typedefs.h"
#include "../task_integrals_common.h"
#include "../../basis/basis.h"

namespace mpqc {
namespace integrals {

/*! \brief Base class for screeners will never skip any integrals.
 *
 * Derived Classes should overload the skip functions which they intend to
 * change.
 */
class Screener {
  public:
    Screener() = default;
    Screener(Screener const &) = default;
    Screener(Screener &&) = default;
    Screener &operator=(Screener &&) = default;
    Screener &operator=(Screener const &) = default;
    virtual ~Screener() = default;

    /// Two loop Outer Screen.
    virtual bool skip(int64_t, Shell const &, ShellVec const &) {
        return false;
    }

    /// Two loop inner screen.
    virtual bool skip(int64_t, int64_t, Shell const &, Shell const &) {
        return false;
    }

    /// Three loop Outer Screen.
    virtual bool
    skip(int64_t, Shell const &, ShellVec const &, ShellVec const &) {
        return false;
    }

    /// Three loop Middle Screen
    virtual bool skip(int64_t, int64_t, Shell const &, Shell const &,
                      ShellVec const &) {
        return false;
    }

    /// Three loop Inner Screen.
    virtual bool skip(int64_t, int64_t, int64_t, Shell const &, Shell const &,
                      Shell const &) {
        return false;
    }

    /// Four loop Outer Screen.
    virtual bool skip(int64_t, Shell const &, ShellVec const &,
                      ShellVec const &, ShellVec const &) {
        return false;
    }

    /// Four loop Second Screen.
    virtual bool skip(int64_t, int64_t, Shell const &, Shell const &,
                      ShellVec const &, ShellVec const &) {
        return false;
    }

    /// Four loop Third Screen.
    virtual bool skip(int64_t, int64_t, int64_t, Shell const &, Shell const &,
                      Shell const &, ShellVec const &) {
        return false;
    }

    /// Four loop Inner Screen.
    virtual bool skip(int64_t, int64_t, int64_t, int64_t, Shell const &,
                      Shell const &, Shell const &, Shell const &) {
        return false;
    }
};

struct init_base_screen {
    template <typename E, unsigned long N>
    Screener operator()(detail::IdxVec const &,
                        detail::ShrBases<N> const &, ShrPool<E> const &) {
        return Screener();
    }
};

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_INTEGRALSCREENERS_H
