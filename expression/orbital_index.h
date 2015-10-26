//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ORBITAL_INDEX_H
#define TILECLUSTERCHEM_ORBITAL_INDEX_H

#include<string>

namespace mpqc{

/* Class to represent Orbital using Char Index
    *   Here's the key index dictionary that can be used (\sa to_space):
    *   MO Orbitals
    *     - m, n -> occ(occupied)
    *     - i,j,k,l -> actocc (actual occupied)
    *     - m', n' -> inactocc (inactive/core occupied)
    *     - x, y -> active (active orbital used in MR)
    *     - a,b,c,d -> virt (virtual)
    *     - p,q,r,s -> any (orbital basis, obs)
    *     - a', b', c', d' -> othervirt (cabs = complete basis - orbital basis)
    *     - A', B', C', D' -> allvirt (complete virtuals = a + a')
    *     - P',Q',R',S' -> allany (complete basis, cbs)
    *
    *   AO Orbitals(greek letter)
    *     - κ λ  μ ν -> obs(orbital basis)
    *     - Κ Λ Μ Ν -> dfbs(density fitting basis)
    *     - α β γ δ -> abs(auxilary basis)
    *
    *     one digit is allowed after letter
    *     for example
    *     i1, i2, a1, a2, P1', a1' ...
 */
class OrbitalIndex{
public:
    enum class Index {
        inactocc = 1,
        actocc = 2,
        occ = 3,
        active = 5,
        virt = 4,
        any = 7,
        othervirt = 6 ,
        allvirt = 10,
        allany = 13,
        obs = -1,
        abs = -2,
        dfbs = -3
    };

    static const wchar_t inactocc_wchar[2];
    static const wchar_t actocc_wchar[2];
    static const wchar_t occ_wchar[2];
    static const wchar_t active_wchar[2];
    static const wchar_t virt_wchar[2];
    static const wchar_t any_wchar[2];
    static const wchar_t othervirt_wchar[2];
    static const wchar_t allvirt_wchar[2];
    static const wchar_t allany_wchar[2];
    static const wchar_t obs_wchar[2];
    static const wchar_t dfbs_wchar[2];
    static const wchar_t abs_wchar[2];

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
