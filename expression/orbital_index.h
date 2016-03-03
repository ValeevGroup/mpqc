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
    *     - i,j,k,l -> actocc (corr occupied)
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
    *     - ρ σ τ υ -> ribs(obs + abs)
    *
    *     digit is allowed after letter
    *     for example
    *     i1, i2, a1, a2, P'1, a'1, i12, a'14 ...
 */
class OrbitalIndex{
public:
    // positive for molecular orbital index
    // negative for atomic orbital index
    enum class Index {
        inactocc = 1,
        active = 2,
        actocc = 3,
        occ = 4,
        virt = 5,
        any = 9,
        othervirt = 10 ,
        allvirt = 15,
        allany = 19,
        obs = -1,
        abs = -2,
        dfbs = -4,
        ribs = -3
    };

    enum class Spin{
        Alpha = 1,
        Beta = -1,
        None = 0
    };

    // constant wchar_t used to map to Index
    static const wchar_t inactocc_wchar[2];
    static const wchar_t actocc_wchar[2];
    static const wchar_t occ_wchar[2];
    static const wchar_t active_wchar[2];
    static const wchar_t virt_wchar[2];
    static const wchar_t any_wchar[2];
    static const wchar_t othervirt_wchar[2];
    static const wchar_t allvirt_wchar[2];
    static const wchar_t allany_wchar[2];
    static const wchar_t obs_wchar[4];
    static const wchar_t dfbs_wchar[4];
    static const wchar_t abs_wchar[4];
    static const wchar_t ribs_wchar[4];

    OrbitalIndex() = default;
    OrbitalIndex(OrbitalIndex const &) = default;
    OrbitalIndex(OrbitalIndex &&) = default;
    OrbitalIndex& operator=(OrbitalIndex const &) = default;
    OrbitalIndex& operator=(OrbitalIndex &&) = default;


    OrbitalIndex(std::wstring letter);
    OrbitalIndex(OrbitalIndex::Index index, OrbitalIndex::Spin spin);

    // if the have the same index
    bool operator==(OrbitalIndex const &) const;
    bool operator==(const OrbitalIndex::Index ) const;
    bool operator<(const OrbitalIndex&)const;
    bool operator>(const OrbitalIndex&)const;

    // if the same index and name
    bool same(const OrbitalIndex& other) const;

    // get index
    const Index &index() const {
        return index_;
    }

    // get spin
    const Spin &spin() const {
        return spin_;
    }

    // get index name
    const std::wstring &name() const {
        return name_;
    }

    // if atomic orbital index
    bool is_ao() const;

    // if molecular orbital index
    bool is_mo() const;

    // if obs
    bool is_mo_in_obs() const;

    // if abs
    bool is_mo_in_abs() const;

    // if ribs
    bool is_mo_in_ribs() const;

    OrbitalIndex mo_to_ao () const;

    std::string to_ta_expression() const;

private:
    void init(const wchar_t *letter);

    Index wchar_to_index(const wchar_t);
    Index wchar_with_prime_to_index(const wchar_t);

private:
    Index index_;
    Spin spin_;
    std::wstring name_;
};

}


#endif //TILECLUSTERCHEM_ORBITAL_INDEX_H
