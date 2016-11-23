#ifndef SRC_MPQC_UTIL_MISC_UNIT_H_
#define SRC_MPQC_UTIL_MISC_UNIT_H_

#include <cmath>
#include <string>

#include "mpqc/util/external/c++/memory"
#include "mpqc/util/misc/exception.h"

namespace mpqc {

namespace detail {
template <typename real_t>
struct FundamentalConstants {
  // fundamental constants proper
  virtual real_t bohr_radius() const = 0;
  virtual real_t elementary_charge() const = 0;
  virtual real_t electron_mass() const = 0;
  virtual real_t Avogadro_constant() const = 0;
  virtual real_t Planck_constant() const = 0;

  // nonfundamental but approved
  virtual real_t atomic_mass_unit() const = 0;

  // derived
  virtual real_t Hartree_to_electron_volt() const = 0;

  // exact in the CODATA system
  static constexpr real_t _electric_constant = 8.854187817e-12;  // F/m
  real_t electric_constant() const {
    return _electric_constant;  // F/m
  }
  // exact
  static constexpr real_t _thermochemical_calorie = 4.184;  // F/m
  real_t thermochemical_calorie() const {
    return _thermochemical_calorie;  // J
  }
  // atomic unit of dipole moment in Debye
  // this is fixed at the value used by MPQC2
  static constexpr real_t _atomic_unit_to_debye = 2.541765;
  real_t atomic_unit_to_debye() const {
    return _atomic_unit_to_debye;  // J
  }
};
}  // namespace detail

namespace constants {

template <typename Real>
struct codata_2014 {
  using real_t = Real;
  static constexpr real_t bohr_radius = 5.2917721067e-11;        // m
  static constexpr real_t elementary_charge = 1.6021766208e-19;  // C
  static constexpr real_t electron_mass = 9.10938356e-31;        // kg
  static constexpr real_t Avogadro_constant = 6.022140857e23;    // mol^-1
  static constexpr real_t Planck_constant = 6.626070040e-34;     // J s
  static constexpr real_t atomic_mass_unit = 1.660539040e-27;    // kg
  // auxiliary conversions
  static constexpr real_t Hartree_to_electron_volt =
      elementary_charge /
      ((4 * M_PI * detail::FundamentalConstants<Real>::_electric_constant) *
       bohr_radius);
};
template <typename Real>
struct codata_2010 {
  using real_t = Real;
  static constexpr real_t bohr_radius = 5.2917721092e-11;       // m
  static constexpr real_t elementary_charge = 1.602176565e-19;  // C
  static constexpr real_t electron_mass = 9.10938291e-31;       // kg
  static constexpr real_t Avogadro_constant = 6.02214129e23;    // mol^-1
  static constexpr real_t Planck_constant = 6.62606957e-34;     // J s
  static constexpr real_t atomic_mass_unit = 1.660538921e-27;   // kg
  // auxiliary conversions
  static constexpr real_t Hartree_to_electron_volt =
      elementary_charge /
      ((4 * M_PI * detail::FundamentalConstants<Real>::_electric_constant) *
       bohr_radius);
};
/// \note based on CODATA1986, except the atomic mass unit
/// see http://physics.nist.gov/cuu/pdf/codata86.pdf
template <typename Real>
struct mpqc2 {
  using real_t = Real;
  static constexpr real_t bohr_radius = 5.29177249e-11;        // m
  static constexpr real_t elementary_charge = 1.60217733e-19;  // C
  static constexpr real_t electron_mass = 9.1093897e-31;       // kg
  static constexpr real_t Avogadro_constant = 6.0221367e23;    // mol^-1
  static constexpr real_t Planck_constant = 6.6260755e-34;     // J s
  // non-CODATA1986
  static constexpr real_t atomic_mass_unit = 1.6605655e-27;  // kg
  // auxiliary conversion explicitly typed in MPQC (not consistent with
  // the above CODATA1986 units)
  static constexpr real_t Hartree_to_electron_volt = 27.2113834;
};

};  // namespace constants

template <typename Data>
struct FundamentalConstants
    : detail::FundamentalConstants<typename Data::real_t> {
  using real_t = typename Data::real_t;
  real_t bohr_radius() const override { return Data::bohr_radius; }
  real_t elementary_charge() const override { return Data::elementary_charge; }
  real_t electron_mass() const override { return Data::electron_mass; }
  real_t Avogadro_constant() const override { return Data::Avogadro_constant; }
  real_t Planck_constant() const override { return Data::Planck_constant; }
  real_t atomic_mass_unit() const override { return Data::atomic_mass_unit; }
  real_t Hartree_to_electron_volt() const override {
    return Data::Hartree_to_electron_volt;
  }
};

class UnitFactory;

/// The Unit class is used to perform unit conversions.
class Unit {
 public:
  ~Unit() = default;

  /// The conversion factor from this to u.
  double to(const Unit& u) const;
  /// The conversion factor from u to this.
  double from(const Unit& u) const;

  /// The conversion factor from this to the corresponding atomic unit.
  double to_atomic_units() const;
  /// The conversion factor from the corresponding atomic unit to this.
  double from_atomic_units() const;

 private:
  friend class UnitFactory;
  friend std::string to_string(const Unit& unit);

  std::string strrep_;
  double to_atomic_units_;
  std::shared_ptr<detail::FundamentalConstants<double>> constants_;

  void parse();

  /// Create using a string representation, like "kcal/mol".
  Unit(const std::string& strrep,
       const std::shared_ptr<detail::FundamentalConstants<double>>& constants);
};

inline std::string to_string(const Unit& unit) { return unit.strrep_; }

class UnitFactory {
 public:
  UnitFactory(std::string version);

  Unit make_unit(const std::string& unit_str) {
    return Unit(unit_str, constants_);
  }

  static std::shared_ptr<UnitFactory> get_default() {
    return instance();
  }
  static void set_default(const std::shared_ptr<UnitFactory>& factory) {
    instance() = factory;
  }

 private:
  std::shared_ptr<detail::FundamentalConstants<double>> constants_;

  static std::shared_ptr<UnitFactory>& instance() {
    static std::shared_ptr<UnitFactory> instance_ = std::make_shared<UnitFactory>("CODATA2014");
    return instance_;
  }

};

}  // namespace mpqc

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
