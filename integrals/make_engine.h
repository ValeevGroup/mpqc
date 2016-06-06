#pragma once
#ifndef MPQC_INTEGRAL_MAKEENGINE_H
#define MPQC_INTEGRAL_MAKEENGINE_H

#include "../basis/basis.h"
// #include "../basis/cluster_shells.h"

#include "task_integrals_common.h"

#include "../molecule/molecule.h"
#include "../molecule/cluster_collapse.h"

#include <mpqc/chemistry/qc/f12/f12_utility.h>

#include <libint2/engine.h>
#include <iostream>

namespace mpqc {
namespace integrals {

template <typename... Bases>
libint2::TwoBodyEngine<libint2::Coulomb> make_2body(Bases &&... basis) {
  int max_am = std::max({basis.max_am()...});
  std::size_t max_nprim = std::max({basis.max_nprim()...});
  return libint2::TwoBodyEngine<libint2::Coulomb>{max_nprim, max_am};
}

template <typename... Bases>
libint2::TwoBodyEngine<libint2::cGTG> make_2body_cGTG(
    std::vector<std::pair<double, double>> const &params, Bases &&... basis) {
  int max_am = std::max({basis.max_am()...});
  std::size_t max_nprim = std::max({basis.max_nprim()...});
  return libint2::TwoBodyEngine<libint2::cGTG>{
      max_nprim, max_am, 0, std::numeric_limits<double>::epsilon(),
      f12::gtg_params_squared(params)};
}

template <typename... Bases>
libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> make_2body_cGTG_C(
    std::vector<std::pair<double, double>> const &params, Bases &&... basis) {
  int max_am = std::max({basis.max_am()...});
  std::size_t max_nprim = std::max({basis.max_nprim()...});
  return libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb>{
      max_nprim, max_am, 0, std::numeric_limits<double>::epsilon(),
      f12::gtg_params_squared(params)};
}

// Function to return the q_vector given a basis
using q_vector = std::vector<std::pair<double, std::array<double, 3>>>;

inline q_vector make_q(molecule::Molecule const &mol) {
  q_vector q;

  for (auto const &cluster : mol) {
    for (auto const &atom : mpqc::molecule::collapse_to_atoms(cluster)) {
      auto const &c = atom.center();
      std::array<double, 3> O = {{c[0], c[1], c[2]}};
      const double charge = atom.charge();

      q.emplace_back(charge, std::move(O));
    }
  }

  return q;
}

inline libint2::Engine
make_1body(std::string const &type, basis::Basis const &bs,
           molecule::Molecule const &mol) {

    // Vector to hold q.
    q_vector q;

    libint2::Operator itype;
    if (type == "overlap") {
        itype = libint2::Operator::overlap;
    } else if (type == "kinetic") {
        itype = libint2::Operator::kinetic;
    } else if (type == "nuclear") {
        itype = libint2::Operator::nuclear;
        q = make_q(mol);
    } else if (type == "emultipole1") {
        itype = libint2::Operator::emultipole1;
    } else {
        std::terminate();
    }

    libint2::Engine engine(itype, bs.max_nprim(),
                                  static_cast<int>(bs.max_am()), 0ul);

    if (itype == libint2::Operator::nuclear) {
        engine.set_params(std::move(q));
    }

    return engine;
}

template <typename... Bases>
inline ShrPool<libint2::TwoBodyEngine<libint2::Coulomb>> make_2body_shr_pool(
    Bases &&... bases) {
  return std::make_shared<Epool<libint2::TwoBodyEngine<libint2::Coulomb>>>(
      make_2body(std::forward<Bases>(bases)...));
}

// template <typename... Bases>
// inline ShrPool<libint2::TwoBodyEngine<libint2::cGTG>> make_2body_cGTG_shr_pool(
//     Bases &&... bases) {
//   // 6 gaussian fit to e^{-0.2 * x} determined in Mathematica
//   std::vector<double> coeffs = {
//       0.68676713614272009761934351651599429330911561677890,
//       0.17924629219445404021936462220200583501183244721279,
//       0.064220682110261371716883231597854661479450500089063,
//       0.033621547173808474311407293716847834106535970827386,
//       0.019592661683693461238951072644210502828610282529384,
//       0.016551662481340943325955850911964947165735196578362};
//   std::vector<double> exps = {
//       0.027451281897186570143120446205078791294303205557980,
//       0.29694066103890138937362272703869530541284211241605,
//       1.3223782109312611325083149412088475466887201014102,
//       4.8167483836252699864760932100830708335301688919489,
//       19.524436889090393389012571590744100348875323294129,
//       159.14994927540994111839660755978500418600010098139};
//   std::vector<std::pair<double, double>> params;
// 
//   for (auto i = 0; i < coeffs.size(); ++i) {
//     params.push_back(std::make_pair(exps[i], coeffs[i]));
//   }
// 
//   return std::make_shared<Epool<libint2::TwoBodyEngine<libint2::cGTG>>>(
//       make_2body_cGTG(std::move(params), std::forward<Bases>(bases)...));
// }

inline ShrPool<libint2::Engine>
make_1body_shr_pool(std::string const &type, basis::Basis const &bs,
                    molecule::Molecule const &mol) {
    return std::make_shared<Epool<libint2::Engine>>(
          make_1body(type, bs, mol));
}

// template <typename... Bases>
// inline ShrPool<libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb>>
// make_2body_cGTG_C_shr_pool(Bases &&... bases) {
//   // 6 gaussian fit to e^{-0.2 * x} determined in Mathematica
//   std::vector<double> coeffs = {
//       0.68676713614272009761934351651599429330911561677890,
//       0.17924629219445404021936462220200583501183244721279,
//       0.064220682110261371716883231597854661479450500089063,
//       0.033621547173808474311407293716847834106535970827386,
//       0.019592661683693461238951072644210502828610282529384,
//       0.016551662481340943325955850911964947165735196578362};
//   std::vector<double> exps = {
//       0.027451281897186570143120446205078791294303205557980,
//       0.29694066103890138937362272703869530541284211241605,
//       1.3223782109312611325083149412088475466887201014102,
//       4.8167483836252699864760932100830708335301688919489,
//       19.524436889090393389012571590744100348875323294129,
//       159.14994927540994111839660755978500418600010098139};
//   std::vector<std::pair<double, double>> params;
// 
//   for (auto i = 0; i < coeffs.size(); ++i) {
//     params.push_back(std::make_pair(exps[i], coeffs[i]));
//   }
// 
//   return std::make_shared<
//       Epool<libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb>>>(
//       make_2body_cGTG_C(std::move(params), std::forward<Bases>(bases)...));
// }

}  // namespace integrals
}  // namespace mpqc

#endif  // MPQC_INTEGRAL_MAKEENGINE_H
