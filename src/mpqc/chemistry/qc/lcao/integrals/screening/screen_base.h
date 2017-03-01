#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCREEN_BASE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCREEN_BASE_H_

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/math/groups/petite_list.h"

#include <tiledarray.h>

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
  virtual bool skip(int64_t);
  virtual bool skip(int64_t) const;

  virtual bool skip(int64_t, int64_t);
  virtual bool skip(int64_t, int64_t) const;

  virtual bool skip(int64_t, int64_t, int64_t);
  virtual bool skip(int64_t, int64_t, int64_t) const;

  virtual bool skip(int64_t, int64_t, int64_t, int64_t);
  virtual bool skip(int64_t, int64_t, int64_t, int64_t) const;

  /*! \brief returns an estimate of shape norms for the given basis vector, in
   * presence of symmetry described by a math::PetiteList object.
   *
   * This function will only compute the estimate for tiles which are
   * considered local by the pmap, the user will be responsible for using the
   * world based constructor for the TA::Shape.
   *
   * \warning The base class will error for all inputs.
   */
  virtual TA::Tensor<float> norm_estimate(
      madness::World &world, std::vector<gaussian::Basis> const &bs_array,
      TA::Pmap const &pmap,
      const math::PetiteList &plist =
          math::SymmPetiteList<math::PetiteList::Symmetry::e>(),
          bool replicate = false) const;
};

}  // namespace  lcao
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_SCREEN_BASE_H_
