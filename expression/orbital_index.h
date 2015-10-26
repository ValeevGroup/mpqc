//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ORBITAL_INDEX_H
#define TILECLUSTERCHEM_ORBITAL_INDEX_H

#include<string>

namespace mpqc{

/* Class to represent Orbital using Char Index
    *   Here's the key index dictionary that can be used (\sa to_space):
    *     - i,j,k,l -> occ (occupied)
    *     - a,b,c,d -> virt (virtual)
    *     - p,q,r,s -> any (orbital basis, obs)
    *     - P',Q',R',S' -> allany (complete basis, cbs)
    *     - a', b', c', d' -> othervirt (cabs = complete basis - orbital basis)
    *     - A', B', C', D' -> allvirt (complete virtuals = a + a')
    *     one digit is allowed after letter
    *     for example
    *     i1, i2, a1, a2, P1', a1' ...
 */
class OrbitalIndex{
public:
    enum class Index {occ = 1, virt = 2, any = 3, othervirt = 4 ,allvirt = 6, allany = 7};

    static const wchar_t occ_char[2];
    static const wchar_t virt_char[2];
    static const wchar_t any_char[2];
    static const wchar_t othervirt_char[2];
    static const wchar_t allvirt_char[2];
    static const wchar_t allany_char[2];

    OrbitalIndex() = default;
    OrbitalIndex(OrbitalIndex const &) = default;
    OrbitalIndex(OrbitalIndex &&) = default;
    OrbitalIndex& operator=(OrbitalIndex const &) = default;
    OrbitalIndex& operator=(OrbitalIndex &&) = default;

    OrbitalIndex(const wchar_t *letter);

    OrbitalIndex(std::wstring letter);

    bool operator==(OrbitalIndex const &);
    bool operator==(const OrbitalIndex::Index );

    const Index &index() const {
        return index_;
    }

private:
    void init(const wchar_t *letter);

private:
    Index index_;
    std::wstring name_;
};


    //TODO frirend function to determine space
}


#endif //TILECLUSTERCHEM_ORBITAL_INDEX_H
