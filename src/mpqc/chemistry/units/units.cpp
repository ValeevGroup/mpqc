
#include "mpqc/chemistry/units/units.h"

#include <cmath>

namespace mpqc {

static inline int eq(const std::unique_ptr<char[]> &a, const char *b) {
  return !strcmp(a.get(), b);
}

Unit::Unit(const std::string &strrep,
           const std::shared_ptr<detail::FundamentalConstants<double>> &consts)
    : strrep_(strrep), constants_(consts) {
  parse();
}

double Unit::to(const Unit &unit) const {
  return to_atomic_units_ / unit.to_atomic_units_;
}

double Unit::from(const Unit &unit) const { return 1.0 / to(unit); }

double Unit::to_atomic_units() const { return to_atomic_units_; }

double Unit::from_atomic_units() const { return 1.0 / to_atomic_units_; }

void Unit::parse() {
  to_atomic_units_ = 1.0;

  // physical constants used for atomic unit conversion factors
  const double a0 = constants_->bohr_radius();        // m
  const double e = constants_->elementary_charge();   // C
  const double me = constants_->electron_mass();      // kg
  const double e0 = constants_->electric_constant();  // F/m
  const double NA = constants_->Avogadro_constant();  // mol-1
  const double h = constants_->Planck_constant();     // J s
  const double hbar = h / (2 * M_PI);                 // J s
  const double c = constants_->speed_of_light();      // m/s

  // derived au conversion factors
  const double Ea = e * e / ((4.0 * M_PI * e0) * a0);  // J

  // other conversions
  const double amu = constants_->atomic_mass_unit();        // kg
  const double cal = constants_->thermochemical_calorie();  // J

  int invert = 0;
  const char *rest = strrep_.c_str();

  while (rest) {
    const char *end = ::strpbrk(rest, " */");
    int nchar;
    if (end) {
      nchar = end - rest;
    } else {
      nchar = strlen(rest);
    }
    std::unique_ptr<char[]> unitstring = std::make_unique<char[]>(nchar + 1);
    memcpy(unitstring.get(), rest, nchar);
    unitstring[nchar] = '\0';

    double factor = 1.0;
    if (eq(unitstring, "bohr") || eq(unitstring, "bohrs")) {
    } else if (eq(unitstring, "hartree") || eq(unitstring, "hartrees") ||
        eq(unitstring, "Hartree") || eq(unitstring, "Hartrees")) {
    } else if (eq(unitstring, "ev") || eq(unitstring, "eV")) {
      factor = 1.0 / constants_->Hartree_to_electron_volt();
    } else if (eq(unitstring, "cm^-1")) {
      factor = (100 * h * c) / Ea;
    } else if (eq(unitstring, "debye")) {
      factor = 1.0 / constants_->atomic_unit_to_debye();
    } else if (eq(unitstring, "radian") || eq(unitstring, "radians")) {
    } else if (eq(unitstring, "mol") || eq(unitstring, "mole")) {
      factor = NA;
    } else if (eq(unitstring, "kcal")) {
      factor = 1000.0 * cal / Ea;
    } else if (eq(unitstring, "kcal_per_mol")) {
      factor = 1000.0 * cal / (Ea * NA);
    } else if (eq(unitstring, "N") || eq(unitstring, "newton")) {
      factor = a0 / Ea;
    } else if (eq(unitstring, "dyne")) {
      factor = 1.0e-5 * a0 / Ea;
    } else if (eq(unitstring, "m") || eq(unitstring, "meter") ||
        eq(unitstring, "meters")) {
      factor = 1.0 / a0;
    } else if (eq(unitstring, "cm") || eq(unitstring, "centimeter")) {
      factor = 1.0e-2 / a0;
    } else if (eq(unitstring, "angstrom") || eq(unitstring, "angstroms") ||
        eq(unitstring, "aangstrom") || eq(unitstring, "aangstroms")) {
      factor = 1.0e-10 / a0;
    } else if (eq(unitstring, "amu")) {
      factor = amu / me;
    } else if (eq(unitstring, "degree") || eq(unitstring, "degrees")) {
      factor = M_PI / 180.0;
    } else if (eq(unitstring, "second") || eq(unitstring, "seconds") ||
        eq(unitstring, "s")) {
      factor = Ea / hbar;
    } else {
      throw InputError("unknown Unit string", __FILE__, __LINE__, "unitstring",
                       unitstring.get());
    }

    unitstring.reset();
    if (invert)
      factor = 1.0 / factor;
    to_atomic_units_ *= factor;
    rest = ::strpbrk(rest, " */");
    while (rest && (*rest == ' ' || *rest == '*' || *rest == '/')) {
      if (*rest == '/')
        invert = !invert;
      rest++;
    }
  }
}

UnitFactory::UnitFactory(std::string system) : system_(std::move(system)) {
  if (system_ == "CODATA2014" || system_ == "2014CODATA") {
    constants_ = std::make_shared<
        FundamentalConstants<constants::codata_2014<double>>>();
  } else if (system_ == "CODATA2010" || system_ == "2010CODATA") {
    constants_ = std::make_shared<
        FundamentalConstants<constants::codata_2010<double>>>();
  } else if (system_ == "CODATA2006" || system_ == "2006CODATA") {
    constants_ = std::make_shared<
        FundamentalConstants<constants::codata_2006<double>>>();
  } else if (system_ == "MPQC2") {
    constants_ =
        std::make_shared<FundamentalConstants<constants::mpqc2<double>>>();
  } else {
    throw mpqc::InputError("Unknown variant of UnitFactory requested", __FILE__,
                           __LINE__, "system", system_.c_str());
  }
}

Unit UnitFactory::make_unit(const std::string &unit_str) const {
  return Unit(unit_str, constants_);
}

std::shared_ptr<const UnitFactory> UnitFactory::get_default() { return instance(); }

void UnitFactory::set_default(const std::string &system) {
  instance() = std::make_shared<const UnitFactory>(system);
}

std::shared_ptr<const UnitFactory> &UnitFactory::instance() {
  static std::shared_ptr<const UnitFactory> instance_ =
      std::make_shared<const UnitFactory>("2014CODATA");
  return instance_;
}

}  // namespace mpqc

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
