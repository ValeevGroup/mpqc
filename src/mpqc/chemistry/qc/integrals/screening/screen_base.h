

#ifndef MPQC_INTEGRALS_INTEGRALSCREENERS_H
#define MPQC_INTEGRALS_INTEGRALSCREENERS_H


#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>
#include <mpqc/chemistry/qc/basis/basis.h>

namespace mpqc {
namespace integrals {

/*! \brief Base class for screeners will never skip any integrals.
 *
 * Derived Classes should overload the skip functions which they intend to
 * change.
 *
 */
class Screener {
  public:
    Screener() = default;
    Screener(Screener const &) = default;
    Screener(Screener &&) = default;
    Screener &operator=(Screener &&) = default;
    Screener &operator=(Screener const &) = default;
    virtual ~Screener() noexcept = default;

    /*! \brief all skips take shell indices and return true or false.
     *
     * The index passed to skip should be the lower bound basis functions for
     * that shell set.  For example if we have shells a, b, c that are all of
     * size 3 then for the integral (a|c) we the indices pasted to skip should
     * be skip(0,5).
     *
     * When deriving this class it is the developers responsiblity to convert
     * function indices into shell indices.  The motivation for this is to
     * simplify the integral kernel code.
     */
    virtual bool skip(int64_t) { return false; }
    virtual bool skip(int64_t) const { return false; }

    virtual bool skip(int64_t, int64_t) { return false; }
    virtual bool skip(int64_t, int64_t) const { return false; }

    virtual bool skip(int64_t, int64_t, int64_t) { return false; }
    virtual bool skip(int64_t, int64_t, int64_t) const { return false; }

    virtual bool skip(int64_t, int64_t, int64_t, int64_t) { return false; }
    virtual bool skip(int64_t, int64_t, int64_t, int64_t) const {
        return false;
    }

};

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_INTEGRALSCREENERS_H
