#ifndef SRC_MPQC_CHEMISTRY_QC_CC_SOLVERS_H_
#define SRC_MPQC_CHEMISTRY_QC_CC_SOLVERS_H_

#include "mpqc/chemistry/qc/cc/tpack.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/core/exception.h"

namespace mpqc {
namespace cc {

/// Solver updates the CC T amplitudes given the current values
/// and the amplitude equation residuals
/// @tparam T the type representing the n-body amplitude set
template <typename T>
class Solver {
 public:
  virtual ~Solver() = default;

  /// Updates amplitudes \c t1 and \c t2 using the residuals \c r1 and \c r2 .
  /// @warning This function assumes that {\c t1 , \c t2 } and {\c r1 , \c r2 }
  ///          are congruent, but does not assume any particular structure.
  ///          Derived classes \em may impose additional assumptions on the
  ///          structure of the arguments.
  /// @param t1 the 1-body amplitude set (current values on input, updated
  /// values on output)
  /// @param t2 the 2-body amplitude set (current values on input, updated
  /// values on output)
  /// @param r1 the 1-body amplitude equation residual set (contents may be
  /// modified)
  /// @param r2 the 2-body amplitude equation residual set (contents may be
  /// modified)
  virtual void update(T& t1, T& t2, const T& r1, const T& r2, double dE = 0) = 0;

  /// Updates amplitudes \c t{1,2,3} using the residuals \c r{1,2,3} .
  /// @warning This function assumes that \c t{1,2,3} and \c r{1,2,3}
  ///          are congruent, but does not assume any particular structure.
  ///          Derived classes \em may impose additional assumptions on the
  ///          structure of the arguments.
  /// @param t1 the 1-body amplitude set (current values on input, updated
  /// values on output)
  /// @param t2 the 2-body amplitude set (current values on input, updated
  /// values on output)
  /// @param t3 the 3-body amplitude set (current values on input, updated
  /// values on output)
  /// @param r1 the 1-body amplitude equation residual set (contents may be
  /// modified)
  /// @param r2 the 2-body amplitude equation residual set (contents may be
  /// modified)
  /// @param r3 the 3-body amplitude equation residual set (contents may be
  /// modified)
  virtual void update(T& t1, T& t2, T& t3, const T& r1, const T& r2, const T& r3) = 0;

  /// Computes the error for the given residuals \c r1 and \c r2 .

  /// The error is defined as the 2-norm per element, i.e.
  /// \f$ \sqrt{||r_1||_2^2 + ||r_2||_2^2}/(\mathrm{size}(r_1) + \mathrm{size}(r_2)) \f$ .
  /// @param[in] r1 the 1-body amplitude equation residual set
  /// @param[in] r2 the 2-body amplitude equation residual set
  /// @return the error
  virtual double error(const T& r1, const T& r2) {
    return std::sqrt((std::pow(norm2(r1),2) + std::pow(norm2(r2),2))) /
        (size(r1) + size(r2));
  }

  /// Tests if the solver has converged
  /// @param target_precision The desired precision for the final CC energy
  /// @param error The value of the error computed from the one- and two-body residuals
  /// @param dE The energy change between two consecutive solver iterations
  /// @return Whether or not the solver has converged
  virtual bool is_converged(double target_precision, double error, double dE) const = 0;
};

/// DIISSolver updates the CC T amplitudes using DIIS
/// @tparam T the type representing the n-body amplitude set
template <typename T>
class DIISSolver : public Solver<T> {
 public:
  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object; it will be queried for the following keywords
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | diis_start | int | 1 | The DIIS extrapolation will begin on the iteration given by this integer. |
   * | n_diis | int | 8 | This specifies maximum number of data sets to retain. |
   * | diis_damp | real | 0.0 | This nonnegative floating point number is used to dampen the DIIS extrapolation. |
   * | diis_ngroup | int | 1 | The number of iterations in a DIIS group. DIIS extrapolation is only used for the first \c diis_group_nstart of these iterations. If \c diis_ngroup is 1 and \c diis_group_nstart is greater than 0, then DIIS will be used on all iterations after and including the start iteration. |
   * | diis_group_nstart | int | 1 | The number of DIIS extrapolations to do at the beginning of an iteration group.  See the documentation for \c diis_ngroup . |
   */
  // clang-format on
  DIISSolver(const KeyVal& kv)
      : diis_(kv.value<int>("diis_start", 1), kv.value<int>("n_diis", 8), kv.value<double>("diis_damp", 0.0),
              kv.value<int>("diis_ngroup", 1), kv.value<int>("diis_group_nstart", 1)) {}
  virtual ~DIISSolver() = default;


  /// Update the amplitudes using update_only() and extrapolate using DIIS.
  /// @warning This function assumes that {\c t1 , \c t2 } and {\c r1 , \c r2 }
  ///          are congruent, but does not assume any particular structure.
  ///          Derived classes \em may impose additional assumptions on the
  ///          structure of the arguments.
  /// @param t1 the 1-body amplitude set (current values on input, updated
  /// values on output)
  /// @param t2 the 2-body amplitude set (current values on input, updated
  /// values on output)
  /// @param r1 the 1-body amplitude equation residual set (contents may be
  /// modified)
  /// @param r2 the 2-body amplitude equation residual set (contents may be
  /// modified)
  void update(T& t1, T& t2, const T& r1, const T& r2, double dE = 0) override {
    update_only(t1, t2, r1, r2, dE);
    T r1_copy = r1;
    T r2_copy = r2;
    TPack<T> r(r1_copy, r2_copy);
    TPack<T> t(t1, t2);
    diis_.extrapolate(t, r);
    t1 = t[0];
    t2 = t[1];
  }

  void update(T& t1, T& t2,  T& t3, const T& r1, const T& r2, const T& r3) override {
    update_only(t1, t2, t3, r1, r2, r3);
    T r1_copy = r1;
    T r2_copy = r2;
    T r3_copy = r3;
    TPack<T> r(r1_copy, r2_copy, r3_copy);
    TPack<T> t(t1, t2, t3);
    diis_.extrapolate(t, r);
    t1 = t[0];
    t2 = t[1];
    t3 = t[2];
  }

  /// Specifies the extra margin for converging energy.
  /// DIIS solver is converged to precision @c eps if the energy between two successive iterations is converged to
  /// @c eps times the margin. Currently it is fixed empirically at 0.1 , i.e. solver is converged if the difference
  /// between two consecutive iterations is less than 1/10th of the target precision.
  /// @return the extra margin for converging energy
  static double precision_margin_energy() {
    return 0.1;
  }

  /// Tests if the solver has converged
  /// @param target_precision The desired precision for the final CC energy
  /// @param error The value of the error computed from the one- and two-body residuals
  /// @param dE The energy change between two consecutive solver iterations
  /// @return Whether or not the solver has converged
  bool is_converged(double target_precision, double error, double dE) const override {
   return (dE <= target_precision * precision_margin_energy() && error <= target_precision);
  }

  /// Resets the DIIS solver; used when switching to a new solver subspace
  void reset() {
    diis_ = decltype(diis_){};
  }

 protected:
  /// this performs the amplitude update only, to be followed up with DIIS
  /// @warning This function assumes that {\c t1 , \c t2 } and {\c r1 , \c r2 }
  ///          are congruent, but does not assume any particular structure.
  ///          Derived classes \em may impose additional assumptions on the
  ///          structure of the arguments.
  virtual void update_only(T& t1, T& t2, const T& r1, const T& r2, double dE=0) {
    throw ProgrammingError("DIISSolver::update_only must be implemented in the derived class",
                           __FILE__, __LINE__);
  }
  virtual void update_only(T& t1, T& t2, T& t3, const T& r1, const T& r2, const T& r3) {
        throw ProgrammingError("DIISSolver::update_only must be implemented in the derived class",
                               __FILE__, __LINE__);
  }


  TA::DIIS<TPack<T>>& diis() { return diis_; }

  private:
  TA::DIIS<TPack<T>> diis_;
 };

}  // namespace cc
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_CC_SOLVERS_H_ */
