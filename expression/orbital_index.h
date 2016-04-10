//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ORBITAL_INDEX_H
#define TILECLUSTERCHEM_ORBITAL_INDEX_H

#include<string>

namespace mpqc{

/**
    \brief Class to represent Orbital using WChar Index
    *
    *   Here's the key index dictionary that can be used
    *
    *  MO Orbitals
    *     - m, n -> occ(occupied)
    *     - i,j,k,l -> corr_occ (correlated occupied)
    *     - m', n' -> frozen_occ (inactive/core occupied)
    *     - x, y -> active (active orbital used in MR)
    *     - a,b,c,d -> virt (virtual)
    *     - p,q,r,s -> any (orbital basis, obs)
    *     - a', b', c', d' -> othervirt (cabs = complete basis - orbital basis)
    *     - A', B', C', D' -> allvirt (complete virtuals = a + a')
    *     - P',Q',R',S' -> allany (complete basis, cbs)
    *
    *   AO Orbitals(greek letter)
    *     - κ λ  μ ν -> obs(primary orbital basis)
    *     - Α Β Γ Δ -> vbs(sencodary orbital basis)
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
    ///
    /// Index types
    /// positive for molecular orbital index
    /// negative for atomic orbital index
    ///
    enum class Index {
        frozen_occ = 1,
        active = 2,
        corr_occ = 3,
        occ = 4,
        virt = 5,
        any = 9,
        othervirt = 10 ,
        allvirt = 15,
        allany = 19,
        obs = -1,
        vbs = -2,
        abs = -3,
        ribs = -4,
        dfbs = -5
    };

    ///
    /// Spin types
    ///
    enum class Spin{
        Alpha = 1,
        Beta = -1,
        None = 0
    };

    ///
    /// constant wchar_t used to map to Index
    ///
    static const wchar_t frozen_occ_wchar[2];
    static const wchar_t corr_occ_wchar[2];
    static const wchar_t occ_wchar[2];
    static const wchar_t active_wchar[2];
    static const wchar_t virt_wchar[2];
    static const wchar_t any_wchar[2];
    static const wchar_t othervirt_wchar[2];
    static const wchar_t allvirt_wchar[2];
    static const wchar_t allany_wchar[2];
    static const wchar_t obs_wchar[4];
    static const wchar_t vbs_wchar[4];
    static const wchar_t dfbs_wchar[4];
    static const wchar_t abs_wchar[4];
    static const wchar_t ribs_wchar[4];

    OrbitalIndex() = default;
    OrbitalIndex(OrbitalIndex const &) = default;
    OrbitalIndex(OrbitalIndex &&) = default;
    OrbitalIndex& operator=(OrbitalIndex const &) = default;
    OrbitalIndex& operator=(OrbitalIndex &&) = default;

    /**
     * Constructor
     * Construct OrbitalIndex wstring
     * @param letter
     * check description of class for mappings
     */

    OrbitalIndex(std::wstring letter);

    /// check equality by comparing index and spin
    bool operator==(OrbitalIndex const &) const;

    /// check equality by index and spin
    bool operator==(const OrbitalIndex::Index, const OrbitalIndex::Spin ) const;

    /// comparison by index and spin
    bool operator<(const OrbitalIndex&)const;

    /// comparison by index and spin
    bool operator>(const OrbitalIndex&)const;

    /// if the same index and name
    bool same(const OrbitalIndex& other) const;

    /// return index
    const Index &index() const {
        return index_;
    }

    /// return spin
    const Spin &spin() const {
        return spin_;
    }

    /// return index name
    const std::wstring &name() const {
        return name_;
    }

    /// if atomic orbital index
    bool is_ao() const;

    /// if molecular orbital index
    bool is_mo() const;

    /// return true if is mo in obs
    bool is_mo_in_obs() const;

    /// return true if is mo in abs
    bool is_mo_in_abs() const;

    /// return true if is mo in ribs
    bool is_mo_in_ribs() const;

    /// Default MO to AO mapping
    /// othervir, allvir, allany -> ribs
    /// everything else -> obs
    OrbitalIndex mo_to_ao () const;

    /// convert name to TiledArray accepted expression
    /// basicly it converts greek letter to english names
    std::string to_ta_expression() const;

private:
    void init(const wchar_t *letter);

    /// convert wchar to index
    Index wchar_to_index(const wchar_t);

    /// convet wchar with prime, for example a', to index
    Index wchar_with_prime_to_index(const wchar_t);

private:
    Index index_;
    Spin spin_;

    /// the name that user passed in from the constructor
    std::wstring name_;
};

}


#endif //TILECLUSTERCHEM_ORBITAL_INDEX_H
