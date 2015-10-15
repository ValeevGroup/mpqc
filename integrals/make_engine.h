#pragma once
#ifndef TCC_INTEGRAL_MAKEENGINE_H
#define TCC_INTEGRAL_MAKEENGINE_H

#include "../include/libint.h"
#include "../basis/basis.h"
// #include "../basis/cluster_shells.h"

#include "../molecule/molecule.h"
#include "../molecule/cluster_collapse.h"

namespace tcc {
namespace integrals {

template <typename... Bases>
libint2::TwoBodyEngine<libint2::Coulomb> make_2body(Bases &&... basis) {
    int max_am = std::max({basis.max_am()...});
    std::size_t max_nprim = std::max({basis.max_nprim()...});
    return libint2::TwoBodyEngine<libint2::Coulomb>{max_nprim, max_am};
}

// Function to return the q_vector given a basis
using q_vector = std::vector<std::pair<double, std::array<double, 3>>>;

// Old basis code
// inline q_vector make_q(mpqc::basis::Basis const &bs) {
//     q_vector q;
// 
//     // Get the groups of clustered shells
//     for (auto const &cluster_shell : bs.cluster_shells()) {
//         // Each group has a reference to its cluster
//         auto const &cluster = cluster_shell.cluster();
//         for (auto const &atom : molecule::collapse_to_atoms(cluster)) {
// 
//             auto const &c = atom.center();
//             std::array<double, 3> O = {{c[0], c[1], c[2]}};
//             const double charge = atom.charge();
// 
//             q.emplace_back(charge, std::move(O));
//         }
//     }
// 
//     return q;
// }

inline q_vector make_q(mpqc::molecule::Molecule const&mol) {
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
make_1body(std::string const &type, mpqc::basis::Basis const &bs, mpqc::molecule::Molecule const &mol) {

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
    } else if (type == "emultipole2") {
        itype = libint2::OneBodyEngine::emultipole2;
    } else {
        std::terminate();
    }

    libint2::OneBodyEngine engine(itype, bs.max_nprim(),
                                  static_cast<int>(bs.max_am()), 0);

    if(itype == libint2::OneBodyEngine::nuclear){
        engine.set_params(std::move(q));
    }

    return engine;
}

} // namespace integrals
} // namespace tcc

#endif // TCC_INTEGRAL_MAKEENGINE_H
