#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_MAKE_ENGINE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_MAKE_ENGINE_H_

#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/integrals/task_integrals_common.h"
#include "mpqc/chemistry/molecule/cluster_collapse.h"
#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/integrals/f12_utility.h"

#include <iostream>
#include <initializer_list>

#include <libint2/engine.h>

namespace mpqc {
namespace lcao {
namespace gaussian {

/// makes an engine for computing integrals of operator \c oper over bases \c
/// bases

///
template <typename Basis, size_t N>
libint2::Engine make_engine(const libint2::Operator &oper,
                            std::array<std::reference_wrapper<Basis>, N> bases,
                            libint2::BraKet braket = libint2::BraKet::invalid,
                            libint2::any oper_params = libint2::any()) {
  // assign default params and braket, if needed
  if (braket == libint2::BraKet::invalid)
    braket = libint2::default_braket(oper);
  if (oper_params.empty()) oper_params = libint2::default_params(oper);

  int max_am = 0;
  size_t max_nprim = 0;
  for (const auto &bs : bases) {
    max_am = std::max(max_am, static_cast<int>(bs.get().max_am()));
    max_nprim = std::max(max_nprim, static_cast<size_t>(bs.get().max_nprim()));
  }
  const auto deriv_order = 0;
  libint2::Engine result{oper, max_nprim, max_am, deriv_order,
                         std::numeric_limits<libint2::real_t>::epsilon(),
                         oper_params, braket};
  return result;
}

// Function to return the q_vector given a basis
using q_vector = std::vector<std::pair<double, std::array<double, 3>>>;

inline q_vector make_q(Molecule const &mol) {
  q_vector q;

  for (auto const &cluster : mol) {
    for (auto const &atom : mpqc::collapse_to_atoms(cluster)) {
      auto const &c = atom.center();
      std::array<double, 3> O = {{c[0], c[1], c[2]}};
      const double charge = atom.charge();

      q.emplace_back(charge, std::move(O));
    }
  }

  return q;
}

template <typename Basis, size_t N>
inline ShrPool<libint2::Engine> make_engine_pool(
    const libint2::Operator &oper,
    std::array<std::reference_wrapper<Basis>, N> bases,
    libint2::BraKet braket = libint2::BraKet::invalid,
    libint2::any oper_params = libint2::any()) {
  // assign default braket, if needed
  if (braket == libint2::BraKet::invalid) {
    braket = libint2::default_braket(oper);
  }
  // assign default params, if needed
  if (oper_params.empty()) {
    oper_params = libint2::default_params(oper);
  }

  return std::make_shared<utility::TSPool<libint2::Engine>>(
      make_engine(oper, bases, braket, oper_params));
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_MAKE_ENGINE_H_
