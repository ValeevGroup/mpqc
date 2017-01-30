

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCREEN_BASE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCREEN_BASE_H_


#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/chemistry/qc/lcao/basis/basis.h"

namespace mpqc {
namespace lcao {

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
  inline virtual bool skip(int64_t) { return false; }
  inline virtual bool skip(int64_t) const { return false; }

  inline virtual bool skip(int64_t, int64_t) { return false; }
  inline virtual bool skip(int64_t, int64_t) const { return false; }

  inline virtual bool skip(int64_t, int64_t, int64_t) { return false; }
  inline virtual bool skip(int64_t, int64_t, int64_t) const { return false; }

  inline virtual bool skip(int64_t, int64_t, int64_t, int64_t) { return false; }
  inline virtual bool skip(int64_t, int64_t, int64_t, int64_t) const {
    return false;
  }

  /*! \brief a debugging function that should validate that the screener is
   * working as expected.
   *
   * \param a function index of first shell
   * \param b function index of second shell
   * \param c function index of third shell
   * \param d function index of fourth shell
   *
   * \param buf is a pointer to the integral buffer
   *
   * \param size should be the number of functions in the buffer
   *
   * \return true if screening was successful
   */
  inline virtual bool validate(int64_t a, int64_t b, int64_t c, int64_t d,
                        double const *buf, int64_t size) const {
    return true;
  }

  /*! \brief a debugging function that should validate that the screener is
   * working as expected.
   *
   *
   * \param a function index of first shell
   * \param b function index of second shell
   * \param c function index of third shell
   *
   * \param buf is a pointer to the integral buffer
   *
   * \param size should be the number of functions in the buffer
   *
   * \return true if screening was successful
   */
  inline virtual bool validate(int64_t a, int64_t b, int64_t c, double const *buf,
                        int64_t size) const {
    return true;
  }
};

}  // namespace  lcao
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCREEN_BASE_H_
