
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_MASSES_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_MASSES_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

namespace mpqc {

/// This singleton provides versioned access to the atomic data.
class AtomicData {
 public:
  struct Isotope {
    std::string label;  // e.g. 12C
    int mass_number;    // mass # = # of protons + # of neutrons
    double mass;        // relative atomic mass (a.m.u.)
    double abundance;   // natural abundance (0..1), (quiet) NAN if not defined
  };

  virtual ~AtomicData() { }

  static std::shared_ptr<const AtomicData> get_default();
  static void set_default(std::string version);

  /// @brief the isotope mass accessor
  ///
  /// @param atomic_number the atomic number (i.e., the number of protons)
  /// @param mass_number the mass number (i.e., the number of protons and
  /// neutrons);
  ///        the default value is "-1", which refers to the most abundant
  ///        isotope.
  /// @return optional value of the relative atomic mass of the specified
  /// isotope,
  ///         in atomic mass units; result is false if the value is not found.
  /// @throws Uncomputable if this atomic number is now known.
  /// @note the isotope masses encoded in isotope_data_ and reported by this
  /// function
  ///       are expressed in a.m.u., i.e. relative to the
  ///       1/12 of the mass of a C atom in its electronic and nuclear ground
  ///       states
  ///       and its rest coordinate system. Therefore they can be converted to
  ///       the atomic
  ///       units safely using the Units library using any system of fundamental
  ///       constants.
  boost::optional<double> isotope_mass(size_t atomic_number,
                                       int64_t mass_number = -1) const;

  /// @return optional value of the mass number of the most abundant isotope;
  /// result is false
  ///         if the most abundant isotope is not known.
  /// @throws Uncomputable if this atomic number is not known.
  boost::optional<const Isotope&> most_abundant_isotope(
      size_t atomic_number) const;

  /// @return a string identifying the version of the atomic data
  virtual std::string version() const = 0;

 protected:
  AtomicData(std::multimap<size_t, Isotope>&& data) : isotope_data_(data) {}

 private:
  const std::multimap<size_t, Isotope> isotope_data_;  // atomic # -> sequence
                                                       // of Isotopes
                                                       // AtomicData() =
                                                       // default;
  static std::shared_ptr<const AtomicData> instance_;

  static std::multimap<size_t, AtomicData::Isotope>::const_iterator
  most_abundant_(
      const std::multimap<size_t, AtomicData::Isotope>::const_iterator& first,
      const std::multimap<size_t, AtomicData::Isotope>::const_iterator& second);
};

/// AtomicData implementation using version 4.1 of <a
/// href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">the
/// NIST database of Atomic Weights and Isotopic Compositions </a>
class NIST_v41_AtomicData : public AtomicData {
 public:
  static const char version_str[];
  NIST_v41_AtomicData() : AtomicData(make_isotope_data()) {}

  std::string version() const override { return version_str; }

 private:
  static std::multimap<size_t, Isotope> make_isotope_data();
};

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_ATOM_MASSES_H_
