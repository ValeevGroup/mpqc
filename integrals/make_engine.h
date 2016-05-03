#pragma once
#ifndef MPQC_INTEGRAL_MAKEENGINE_H
#define MPQC_INTEGRAL_MAKEENGINE_H

#include "../basis/basis.h"
// #include "../basis/cluster_shells.h"

#include "task_integrals_common.h"

#include "../molecule/molecule.h"
#include "../molecule/cluster_collapse.h"

#include <libint2/engine.h>

namespace mpqc {
namespace integrals {

template <typename... Bases>
libint2::TwoBodyEngine<libint2::Coulomb> make_2body(Bases &&... basis) {
    int max_am = std::max({basis.max_am()...});
    std::size_t max_nprim = std::max({basis.max_nprim()...});
    return libint2::TwoBodyEngine<libint2::Coulomb>{max_nprim, max_am};
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

// q must be set a later time since I am not sure about the ordering aspect.
inline libint2::OneBodyEngine
make_1body(std::string const &type, basis::Basis const &bs,
           molecule::Molecule const &mol) {

    // Vector to hold q.
    q_vector q;

    libint2::OneBodyEngine::operator_type itype;
    if (type == "overlap") {
        itype = libint2::OneBodyEngine::overlap;
    } else if (type == "kinetic") {
        itype = libint2::OneBodyEngine::kinetic;
    } else if (type == "nuclear") {
        itype = libint2::OneBodyEngine::nuclear;
        q = make_q(mol);
    } else if (type == "emultipole1") {
        itype = libint2::OneBodyEngine::emultipole1;
    } else {
        std::terminate();
    }

    libint2::OneBodyEngine engine(itype, bs.max_nprim(),
                                  static_cast<int>(bs.max_am()), 0ul);

    if (itype == libint2::OneBodyEngine::nuclear) {
        engine.set_params(std::move(q));
    }

    return engine;
}

inline ShrPool<libint2::OneBodyEngine>
make_1body_shr_pool(std::string const &type, basis::Basis const &bs,
                    molecule::Molecule const &mol) {
    return std::make_shared<Epool<libint2::OneBodyEngine>>(
          make_1body(type, bs, mol));
}

template <typename... Bases>
inline ShrPool<libint2::TwoBodyEngine<libint2::Coulomb>>
make_2body_shr_pool(Bases &&... bases) {
    return std::make_shared<Epool<libint2::TwoBodyEngine<libint2::Coulomb>>>(
          make_2body(std::forward<Bases>(bases)...));
}

} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRAL_MAKEENGINE_H
