#include "mpqc/chemistry/molecule/atomic_data.h"

#include <cmath>

#include "mpqc/util/misc/exception.h"

namespace mpqc {

std::shared_ptr<const AtomicData> AtomicData::instance_ = nullptr;

std::shared_ptr<const AtomicData> AtomicData::get_default() {
  if (!instance_) instance_ = std::make_shared<const NIST_v41_AtomicData>();
  return instance_;
}

void AtomicData::set_default(std::string version) {
  if (version != instance_->version()) {
    if (version == NIST_v41_AtomicData::version_str)
      instance_ = std::make_shared<const NIST_v41_AtomicData>();
  }
}

boost::optional<double> AtomicData::isotope_mass(size_t atomic_number,
                                                 int64_t mass_number) const {
  auto range = isotope_data_.equal_range(atomic_number);
  if (range.second == range.first)
    throw Uncomputable("did not find isotopes for this atomic number", __FILE__,
                       __LINE__);

  auto the_iter =
      (mass_number == -1)
          ? most_abundant_(range.first, range.second)
          : std::find_if(range.first, range.second,
                         [=](const std::pair<size_t, Isotope>& a) -> bool {
                           assert(a.first == atomic_number);
                           return a.second.mass_number == mass_number;
                         });

  typedef boost::optional<double> result_t;
  return (the_iter == range.second) ? result_t()
                                    : result_t(the_iter->second.mass);
}

boost::optional<const AtomicData::Isotope&> AtomicData::most_abundant_isotope(
    size_t atomic_number) const {
  auto range = isotope_data_.equal_range(atomic_number);
  if (range.second == range.first)
    throw Uncomputable("did not find isotopes for this atomic number", __FILE__,
                       __LINE__);
  auto the_iter = most_abundant_(range.first, range.second);

  // return empty optional if no isotope has known abundance
  typedef boost::optional<const AtomicData::Isotope&> result_t;
  return std::isnan(the_iter->second.abundance) ? result_t()
                                                : result_t(the_iter->second);
}

std::multimap<size_t, AtomicData::Isotope>::const_iterator
AtomicData::most_abundant_(
    const std::multimap<size_t, AtomicData::Isotope>::const_iterator& first,
    const std::multimap<size_t, AtomicData::Isotope>::const_iterator& second) {
  return std::max_element(
      first, second, [](const std::pair<size_t, Isotope>& a,
                        const std::pair<size_t, Isotope>& b) -> bool {
        const auto& a_ab = a.second.abundance;
        const auto& b_ab = b.second.abundance;
        if (std::isnan(
                a_ab)) {  // neither has abundance or only b has it? b wins
          return true;
        } else {
          if (std::isnan(b_ab))  // only a has abundance? a wins
            return false;
          else
            return a_ab < b_ab;
        }
      });
}

const char NIST_v41_AtomicData::version_str[] = "NIST-4.1";

std::multimap<size_t, AtomicData::Isotope>
NIST_v41_AtomicData::make_isotope_data() {
  return {
      // clang-format off
      //{Z, {symbol,  A,      mass         , abundance
      {1,   {"H",       1,      1.00782503223, 0.999885}},
      {1,   {"D",       2,      2.01410177812, 0.000115}},
      {1,   {"T",       3,      2.01410177812, NAN}},
      {2,   {"3He",     3,      3.0160293201,  0.00000134}},
      {2,   {"4He",     4,      4.00260325413, 0.00000134}},
      {3,   {"6Li",     6,      6.0151228874,  0.0759}},
      {3,   {"7Li",     7,      7.0160034366,  0.9241}},
      {4,   {"9Be",     9,      9.012183065,   1}},
      {5,   {"10B",    10,     10.01293695,    0.199}},
      {5,   {"11B",    11,     11.00930536,    0.801}},
      {6,   {"12C",    12,     12.0000000,     0.9893}},
      {6,   {"13C",    13,     13.00335483507, 0.0107}},
      {6,   {"14C",    14,     14.0032419884,  NAN}},
      {7,   {"14N",    14,     14.00307400443, 0.99636}},
      {7,   {"15N",    15,     15.00010889888, 0.00364}},
      {8,   {"16O",    16,     15.99491461957, 0.99757}},
      {8,   {"17O",    17,     16.99913175650, 0.00038}},
      {8,   {"18O",    18,     17.99915961286, 0.00205}},
      {9,   {"19F",    19,     18.99840316273, 1}},
      {10,  {"20Ne",   20,     19.9924401762,  0.9048}},
      {10,  {"21Ne",   21,     20.993846685,   0.0027}},
      {10,  {"22Ne",   22,     21.991385114,   0.0925}}
      // clang-format on
  };
}

}  // namespace mpqc
