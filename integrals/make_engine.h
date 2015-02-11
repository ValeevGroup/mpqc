#pragma once
#ifndef TCC_INTEGRAL_MAKEENGINE_H
#define TCC_INTEGRAL_MAKEENGINE_H

#include "../include/libint.h"
#include "../basis/basis.h"

namespace tcc {
namespace integrals {

template <typename... Bases>
libint2::TwoBodyEngine<libint2::Coulomb> make_2body(Bases &&... basis) {
    int max_am = std::max({basis.max_am()...});
    auto max_nprim = std::max({basis.max_nprim()...});
    return libint2::TwoBodyEngine<libint2::Coulomb>{max_nprim, max_am};
}

// q must be set a later time since I am not sure about the ordering aspect.
libint2::OneBodyEngine
make_1body(std::string const &type, basis::Basis const &bs) {

    libint2::OneBodyEngine::integral_type itype;
    if(type == "overlap"){
        itype = libint2::OneBodyEngine::overlap;
    }
    else if(type == "nuclear"){
        itype = libint2::OneBodyEngine::nuclear;
    }
    else if(type == "kinetic"){
        itype = libint2::OneBodyEngine::kinetic;
    } else {
        std::terminate();
    }

    return libint2::OneBodyEngine{itype, bs.max_nprim(),
                                  static_cast<int>(bs.max_am()), 0};
}

} // namespace integrals
} // namespace tcc

#endif // TCC_INTEGRAL_MAKEENGINE_H
