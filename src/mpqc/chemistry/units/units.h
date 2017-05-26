#ifndef MPQC4_SRC_MPQC_CHEMISTRY_UNITS_UNITS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_UNITS_UNITS_H_

#include <cmath>
#include <memory>
#include <string>

#include "mpqc/util/misc/exception.h"

namespace mpqc {

namespace detail {
/// Abstract base for a fundamental constants system
/// @tparam the real type used to report the values
template <typename real_t>
struct FundamentalConstants {
  /// @return a string describing the particular fundamental constants system
  virtual const char* description() const = 0;

  ////////////////////////////////////////////////////////////////////////////
  // fundamental constants proper

  /// @return the value of Bohr radius (m)
  virtual real_t bohr_radius() const = 0;
  /// @return the value of elementary charge (C)
  virtual real_t elementary_charge() const = 0;
  /// @return the value of electron mass (kg)
  virtual real_t electron_mass() const = 0;
  /// @return the value of Avogadro number ( mol\f$^{-1}\f$ )
  virtual real_t Avogadro_constant() const = 0;
  /// @return the value of Planck constant ( J s )
  virtual real_t Planck_constant() const = 0;

  /// @return the atomic mass unit, a.m.u. (kg)
  /// @note nonfundamental but approved in CODATA
  virtual real_t atomic_mass_unit() const = 0;

  ////////////////////////////////////////////////////////////////////////////
  // derived

  /// @return the value of eV / Hartree ratio (unitless)
  virtual real_t Hartree_to_electron_volt() const = 0;

  ////////////////////////////////////////////////////////////////////////////
  // exact in the CODATA system

  static constexpr real_t _electric_constant = 8.854187817e-12;  // F/m
  /// @return the value of electric constant (F / m)
  /// @note this constant is exact in the SI system
  real_t electric_constant() const {
    return _electric_constant;  // F/m
  }
  static constexpr real_t _speed_of_light = 299792458;  // m/s
  /// @return the value of the speed of light in vacuum (m/s)
  /// @note this constant is exact in the SI system
  real_t speed_of_light() const {
    return _speed_of_light;  // m/s
  }
  /// @return the value of thermochemical calorie (J)
  /// @note this is an exact conversion factor
  static constexpr real_t _thermochemical_calorie = 4.184;  // J
  real_t thermochemical_calorie() const {
    return _thermochemical_calorie;  // J
  }
  /// @return the atomic unit of dipole moment (Debye)
  /// @note this is fixed at the value used by MPQC2
  static constexpr real_t _atomic_unit_to_debye = 2.541765;
  real_t atomic_unit_to_debye() const {
    return _atomic_unit_to_debye;  // J
  }
};
}  // namespace detail

namespace constants {

/// The 2014 CODATA revision
template <typename Real>
struct codata_2014 {
  static constexpr const char* const description =
      "2014 CODATA revision: DOI 10.1103/RevModPhys.88.035009";
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

/// The 2010 CODATA revision
template <typename Real>
struct codata_2010 {
  static constexpr const char* const description =
      "2010 CODATA revision: DOI 10.1103/RevModPhys.84.1527";
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
/// The 2006 CODATA revision
/// \note 2006 CODATA set is the default in Gaussian09 (see http://www.gaussian.com/g_tech/g_ur/k_constants.htm)
template <typename Real>
struct codata_2006 {
  static constexpr const char* const description =
      "2006 CODATA revision: DOI 10.1103/RevModPhys.80.633";
  using real_t = Real;
  static constexpr real_t bohr_radius = 5.2917720859e-11;       // m
  static constexpr real_t elementary_charge = 1.602176487e-19;  // C
  static constexpr real_t electron_mass = 9.10938215e-31;       // kg
  static constexpr real_t Avogadro_constant = 6.02214179e23;    // mol^-1
  static constexpr real_t Planck_constant = 6.62606896e-34;     // J s
  static constexpr real_t atomic_mass_unit = 1.660538782e-27;   // kg
  // auxiliary conversions
  static constexpr real_t Hartree_to_electron_volt =
      elementary_charge /
      ((4 * M_PI * detail::FundamentalConstants<Real>::_electric_constant) *
       bohr_radius);
};
/// The constants system used by MPQC2
/// \note based on CODATA1986, except the atomic mass unit
/// see http://physics.nist.gov/cuu/pdf/codata86.pdf
template <typename Real>
struct mpqc2 {
  static constexpr const char* const description = "MPQC 2.3 constants, based on 1986 CODATA";
  using real_t = Real;
  static constexpr real_t bohr_radius = 5.29177249e-11;        // m
  static constexpr real_t elementary_charge = 1.60217733e-19;  // C
  static constexpr real_t electron_mass = 9.1093897e-31;       // kg
  static constexpr real_t Avogadro_constant = 6.0221367e23;    // mol^-1
  static constexpr real_t Planck_constant = 6.6260755e-34;     // J s
  // non-CODATA1986
  static constexpr real_t atomic_mass_unit = 1.6605655e-27;    // kg
  // auxiliary conversion explicitly typed in MPQC (not consistent with
  // the above CODATA1986 units)
  static constexpr real_t Hartree_to_electron_volt = 27.2113834;
};

};  // namespace constants

/// The set of fundamental constants described by \c System
/// @tparam System a data type describing the system of fundamental constants (see namespace ::mpqc::constants )
template <typename System>
struct FundamentalConstants
    : detail::FundamentalConstants<typename System::real_t> {
  const char* description() const override { return System::description; }
  using real_t = typename System::real_t;
  real_t bohr_radius() const override { return System::bohr_radius; }
  real_t elementary_charge() const override { return System::elementary_charge; }
  real_t electron_mass() const override { return System::electron_mass; }
  real_t Avogadro_constant() const override { return System::Avogadro_constant; }
  real_t Planck_constant() const override { return System::Planck_constant; }
  real_t atomic_mass_unit() const override { return System::atomic_mass_unit; }
  real_t Hartree_to_electron_volt() const override {
    return System::Hartree_to_electron_volt;
  }
};

class UnitFactory;

/// \brief The Unit class is used to perform unit conversions.
/// \note Unit conversion factors depend on the fundamental physical constants,
/// hence each Unit object refers to a particular constants system.
/// To ensure consistent use of Unit objects they can only be created using
/// an object of UnitFactory class.
/// \sa UnitFactory
class Unit {
 public:
  ~Unit() { }

  /// The conversion factor from this to \c u .
  double to(const Unit& u) const;
  /// The conversion factor from \c u to this.
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

/// UnitFactory produces Unit objects that refer to a particular system of
/// fundamental constants.

/// UnitFactory is a helper class to ensure that all unit conversions are consistent,
/// i.e. refer to the same system of fundamental constants. Since fundamental constants
/// are updated every few years, it is mandatory to use UnitFactory to produce Unit objects,
/// rather than create them directly. Since typically only 1 UnitFactory needs to be used
/// in the entire program, the UnitFactory should be used as a singleton:
/// \code
/// int main() {
///   UnitFactory::set_default("CODATA2010");  // set the 2010 CODATA revision as default
///   ...
///   auto angstrom = UnitFactory::get_default().make_unit("angstrom");
///
/// }
/// \endcode
class UnitFactory {
 public:
  /// Creates a UnitFactory object corresponding to the given fundamental constants
  /// system.
  /// \param system specifies the fundamental constants system, the allowed values are:
  ///    - "2014CODATA" : the 2014 revision of the fundamental constants (see DOI 10.1103/RevModPhys.88.035009 )
  ///    - "2010CODATA" : the 2010 revision of the fundamental constants (see DOI 10.1103/RevModPhys.84.1527 )
  ///    - "2006CODATA" : the 2006 revision of the fundamental constants (see DOI 10.1103/RevModPhys.80.633 ) (the default of Gaussian09)
  ///    - "MPQC2" : the values of constants used by MPQC version 2.3
  ///
  /// The default is currently "2014CODATA", and may be revised in the future.
  UnitFactory(std::string system);

  ~UnitFactory() { }

  /// \return the name of the fundamental constants system
  const std::string& system() const { return system_; }

  /// \return a shared_ptr to the fundamental constants system object
  std::shared_ptr<const detail::FundamentalConstants<double>> constants() const { return constants_; }

  /// \return the description of the fundamental constants system
  const char* description() const { return constants_->description(); }

  /// makes a Unit object from a given string specification, using the system
  /// of fundamental constants specified by this factory.
  /// \param unit_str the string specification of the desired Unit. See \ref units section.
  /// \return the Unit object
  Unit make_unit(const std::string& unit_str) const;

  /// \return the singleton UnitFactory object (default initialized with 2014CODATA fundamental constants)
  static std::shared_ptr<const UnitFactory> get_default();
  /// sets the singleton UnitFactory object
  /// \param system specifies the fundamental constants system, must be
  /// acceptable as the input to UnitFactory::UnitFactory() ctor
  static void set_default(const std::string& system);

 private:
  std::string system_;
  std::shared_ptr<detail::FundamentalConstants<double>> constants_;

  static std::shared_ptr<const UnitFactory>& instance();

};

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_UNITS_UNITS_H_

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
