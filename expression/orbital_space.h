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
    OrbitalSpace(const OrbitalIndex& index, const Array& tarray)
            : index_(index), coefs_(tarray)
    {}

    ~OrbitalSpace()= default;

    OrbitalIndex &key(){
        return index_;
    }

    const OrbitalIndex &key() const {
        return index_;
    }

    // interface to TA::Array () function
    TA::expressions::TsrExpr<Array,true>
            operator()(const std::string& vars){
        return coefs_(vars);
    };

    // interface to TA::Array () function
    TA::expressions::TsrExpr<const Array,true>
    operator()(const std::string& vars) const {
        return coefs_(vars);
    };

private:

    OrbitalIndex index_;
    Array coefs_;

};

}

#endif //TILECLUSTERCHEM_ORBITAL_SPACE_H
