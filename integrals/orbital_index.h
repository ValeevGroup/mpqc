//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ORBITAL_INDEX_H
#define TILECLUSTERCHEM_ORBITAL_INDEX_H

#include <string>

namespace mpqc{

/* Class to represent Orbital using Char Index
    *   Here's the key index dictionary that can be used (\sa to_space):
    *     - i,j,k,l -> occ (occupied)
    *     - a,b,c,d -> virt (virtual)
    *     - p,q,r,s -> any (orbital basis, obs)
    *     - P',Q',R',S' -> allany (complete basis, cbs)
    *     - a', b', c', d' -> othervirt (cabs = complete basis - orbital basis)
    *     - A', B', C', D' -> allvirt (complete virtuals = a + a')
 */
class OrbitalIndex{
public:
    enum class Index {occ = 1, virt = 2, any = 3, othervirt = 4 ,allvirt = 6, allany = 7};

    static const char occ_char[2];
    static const char virt_char[2];
    static const char any_char[2];
    static const char othervirt_char[2];
    static const char allvirt_char[2];
    static const char allany_char[2];

    OrbitalIndex() = default;
    OrbitalIndex(OrbitalIndex const &) = default;
    OrbitalIndex(OrbitalIndex &&) = default;
    OrbitalIndex& operator=(OrbitalIndex const &) = default;
    OrbitalIndex& operator=(OrbitalIndex &&) = default;

    OrbitalIndex(const char *letter);

    OrbitalIndex(const std::string letter) : OrbitalIndex(letter.c_str()) {}

    bool operator==(OrbitalIndex const &);

private:

    Index index_;
};

    // set the range of key index to use
    const char OrbitalIndex::occ_char[2] = {'i','l'};
    const char OrbitalIndex::virt_char[2] = {'a','d'};
    const char OrbitalIndex::any_char[2] = {'p','s'};
    const char OrbitalIndex::othervirt_char[2] = {'a','d'};
    const char OrbitalIndex::allvirt_char[2] = {'A','D'};
    const char OrbitalIndex::allany_char[2] = {'P','S'};

}


#endif //TILECLUSTERCHEM_ORBITAL_INDEX_H
