//
// Created by Chong Peng on 2/16/16.
//

#ifndef TILECLUSTERCHEM_ORBITAL_SPACE_H
#define TILECLUSTERCHEM_ORBITAL_SPACE_H

#include <memory>

#include "orbital_index.h"
#include "operation.h"
#include "../basis/basis.h"

namespace mpqc{

template <typename Array>
class OrbitalSpace{
public:

    OrbitalSpace() = default;
    OrbitalSpace(const OrbitalIndex& index, const basis::Basis& basis_set, const Array& tarray)
            : mo_index_(index), basis_set_(std::make_shared<basis::Basis>(basis_set)), coefs_(tarray)
    {}

    ~OrbitalSpace()= default;

    OrbitalIndex key(){
        return mo_index_;
    }

private:

    OrbitalIndex mo_index_;
    std::shared_ptr<basis::Basis> basis_set_;
    Array coefs_;

};




}

#endif //TILECLUSTERCHEM_ORBITAL_SPACE_H
