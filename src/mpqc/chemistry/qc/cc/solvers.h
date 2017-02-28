#ifndef SRC_MPQC_CHEMISTRY_QC_CC_SOLVERS_H_
#define SRC_MPQC_CHEMISTRY_QC_CC_SOLVERS_H_

#include "diis.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace cc {

/// Solver updates the CC T amplitudes given the current values
/// and the amplitude equation residuals
/// @tparam T1 the type representing the 1-body amplitude set
/// @tparam T2 the type representing the 2-body amplitude set
template <typename T1, typename T2>
class Solver {
 public:
  virtual ~Solver() = default;

  /// @param t1 the 1-body amplitude set (current values on input, updated
  /// values on output)
  /// @param t2 the 2-body amplitude set (current values on input, updated
  /// values on output)
  /// @param r1 the 1-body amplitude equation residual set (contents may be
  /// modified)
  /// @param r2 the 2-body amplitude equation residual set (contents may be
  /// modified)
  virtual void update(T1& t1, T2& t2, const T1& r1, const T2& r2) = 0;
};

/// DIISSolver updates the CC T amplitudes using DIIS
/// @tparam T1 the type representing the 1-body amplitude set
/// @tparam T2 the type representing the 2-body amplitude set
template <typename T1, typename T2>
class DIISSolver : public Solver<T1, T2> {
 public:
  DIISSolver(const KeyVal& kv)
      : diis_(kv.value<int>("diis_strt", 1), kv.value<int>("n_diis", 8), 0.0,
              kv.value<int>("diis_ngr", 2), kv.value<int>("ngrdiis", 1)) {}
  virtual ~DIISSolver() = default;

  /// Update the amplitudes using update_only() and extrapolate using DIIS.
  /// @param t1 the 1-body amplitude set (current values on input, updated
  /// values on output)
  /// @param t2 the 2-body amplitude set (current values on input, updated
  /// values on output)
  /// @param r1 the 1-body amplitude equation residual set (contents may be
  /// modified)
  /// @param r2 the 2-body amplitude equation residual set (contents may be
  /// modified)
  void update(T1& t1, T2& t2, const T1& r1, const T2& r2) override {
    update_only(t1, t2, r1, r2);
    T1 r1_copy = r1;
    T1 r2_copy = r2;
    T1T2<T1, T2> r(r1_copy, r2_copy);
    T1T2<T1, T2> t(t1, t2);
    diis_.extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
    std::cout << "t1 (after DIIS)" << t1 << std::endl;
    std::cout << "t2 (after DIIS)" << t2 << std::endl;
  }

 protected:
  /// this performs the amplitude update only, without
  virtual void update_only(T1& t1, T2& t2, const T1& r1, const T2& r2) = 0;

 private:
  TA::DIIS<T1T2<T1, T2>> diis_;
};

}  // namespace cc
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_CC_SOLVERS_H_ */
