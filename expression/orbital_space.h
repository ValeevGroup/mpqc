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
    OrbitalSpace(const OrbitalIndex& mo_index, const OrbitalIndex& ao_index, const Array& tarray)
            : mo_index_(mo_index), ao_index_(ao_index), coefs_(tarray)
    {}

    ~OrbitalSpace()= default;

    OrbitalIndex &mo_key(){
        return mo_index_;
    }

    const OrbitalIndex &mo_key() const {
        return mo_index_;
    }

    OrbitalIndex &ao_key(){
        return ao_index_;
    }

    const OrbitalIndex &ao_key() const {
        return ao_index_;
    }

    Array& array() {
        return coefs_;
    }

    const Array& array() const {
        return coefs_;
    }

    /// interface to TA::Array () function
    TA::expressions::TsrExpr<Array, true>
            operator()(const std::string& vars){
        return coefs_(vars);
    };

    /// interface to TA::Array () function
    TA::expressions::TsrExpr<const Array,true>
    operator()(const std::string& vars) const {
        return coefs_(vars);
    };

private:

    OrbitalIndex mo_index_;
    OrbitalIndex ao_index_;
    Array coefs_;

};

}

#endif //TILECLUSTERCHEM_ORBITAL_SPACE_H
