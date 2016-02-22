//
// Created by Chong Peng on 2/16/16.
//

#ifndef TILECLUSTERCHEM_ORBITAL_SPACE_REGISTRY_H
#define TILECLUSTERCHEM_ORBITAL_SPACE_REGISTRY_H

#include "formula_registry.h"
#include "orbital_space.h"


namespace mpqc{

template <typename Array>
class OrbitalSpaceRegistry : public Registry<OrbitalIndex,OrbitalSpace<Array>>{
public:

    using Key = OrbitalIndex;
    using Value = OrbitalSpace<Array>;
    using value_type = typename Registry<Key,Value>::value_type;
    using element_type = typename Registry<Key,Value>::element_type;
    using iterator = typename Registry<Key,Value>::iterator;
    using const_iterator = typename Registry<Key,Value>::const_iterator;

    OrbitalSpaceRegistry()= default;
    OrbitalSpaceRegistry(const element_type& map) : Registry<Key,Value>(map){}
};

} // end of namespace mpqc


#endif //TILECLUSTERCHEM_ORBITAL_SPACE_REGISTRY_H
