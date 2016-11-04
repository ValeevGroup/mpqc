//
// Created by Chong Peng on 10/14/15.
//

#ifndef MPQC_ORBITAL_INDEX_H
#define MPQC_ORBITAL_INDEX_H

#include <string>

#include <TiledArray/error.h>

#include "mpqc/util/misc/string.h"

namespace mpqc {

/**
    \brief OrbitalIndex refers to an OrbitalSpace.
    *
    *   OrbitalIndex is a label (wide-char string) whose semantics follow
    *   a convention common in the molecular electronic structure.
    *   Thus it forms a basis for a language describing mathematical expressions
    *   encountered in this domain.
    *
    *   The labels consist of a single wchar, followed by an optional prime,
    *   followed by zero or more digits, followed by an optional spin label.
    *   The non-digit characters must belong to
    *   the following dictionary:

        <table border="1">

        <tr> <td> Label base <td> OrbitalIndex::Type <td> Description

        <tr><td><tt>m,n</tt><td><tt>occ</tt><td> (active and core) occupied </tr>
        <tr><td><tt>m', n'</tt><td><tt>frozen_occ</tt><td> core occupied </tr>
        <tr><td><tt>i,j,k,l</tt><td><tt>corr_occ</tt><td> active occupied </tr>
        <tr><td><tt>x, y</tt><td><tt>active</tt><td> active (e.g. active space in MR theories) </tr>
        <tr><td><tt>a,b,c,d</tt><td><tt>virt</tt><td> unoccupied/virtual </tr>
        <tr><td><tt>p,q,r,s</tt><td><tt>any</tt><td> occupied or unoccupied expressed in the orbital AO basis and/or virtual AO basis </tr>
        <tr><td><tt>a', b', c', d'</tt><td><tt>othervirt</tt><td> unoccupied orbitals orthogonal to <tt>any</tt> (in F12 theory denoted as CABS) </tr>
        <tr><td><tt>A', B', C', D'</tt><td><tt>allvirt</tt><td> all unoccupied orbitals, i.e. union of <tt>virt</tt> and <tt>othervirt</tt> </tr>
        <tr><td><tt>P', Q', R', S'</tt><td><tt>allany</tt><td> union of <tt>any</tt> and <tt>othervirt</tt> </tr>
        <tr><td><tt>κ, λ, μ, ν</tt><td><tt>obs</tt><td> orbital AO basis, supports at least occupied orbitals, or more usually supports <tt>any</tt> and all of its subsets </tr>
        <tr><td><tt>Α, Β, Γ, Δ</tt><td><tt>vbs</tt><td> unoccupied/virtual AO basis, supports <tt>virt</tt> in <em>dual-basis</em> theories </tr>
        <tr><td><tt>Κ, Λ, Μ, Ν</tt><td><tt>dfbs</tt><td> density-fitting AO basis </tr>
        <tr><td><tt>α, β, γ, δ</tt><td><tt>abs</tt><td> auxiliary AO basis used in F12 theories </tr>
        <tr><td><tt>ρ, σ, τ, υ</tt><td><tt>ribs</tt><td> union of <tt>obs</tt> and <tt>abs</tt>, supports <tt>othervirt</tt>, <tt>allvirt</tt>, and <tt>allany</tt> </tr>

        </table>
 */
class OrbitalIndex {
 public:
  ///
  /// Index types
  /// positive for molecular orbital index
  /// negative for atomic orbital index
  ///
  enum class Type {
    frozen_occ = 1,
    active = 2,
    corr_occ = 3,
    occ = 4,
    virt = 5,
    any = 9,
    othervirt = 10,
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
  enum class Spin { Alpha = 1, Beta = -1, None = 0 };

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
  OrbitalIndex &operator=(OrbitalIndex const &) = default;
  OrbitalIndex &operator=(OrbitalIndex &&) = default;

  /**
   * \brief constructs from a label
   *
   * Construct OrbitalIndex from a std::wstring
   * @param symbol
   * check description of class for mappings
   */
  template <typename String,
            typename = typename std::enable_if<
                not std::is_same<typename std::decay<String>::type, OrbitalIndex>::value>::type>
  OrbitalIndex(String &&symbol);

  /// check equality by comparing index and spin
  bool operator==(OrbitalIndex const &) const;

  /// check inequality by comparing index and spin
  bool operator!=(OrbitalIndex const &) const;

  /// comparison by index and spin
  bool operator<(const OrbitalIndex &) const;

  /// comparison by index and spin
  bool operator>(const OrbitalIndex &) const;

  /// if the same index and name
  bool same(const OrbitalIndex &other) const;

  /// return index
  const Type &index() const { return index_; }

  /// return spin
  const Spin &spin() const { return spin_; }

  /// return index name
  const std::wstring &name() const { return name_; }

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
  OrbitalIndex mo_to_ao() const;

  /// convert name to TiledArray accepted expression,
  /// i.e. it converts greek letters with english equivalents
  std::string to_ta_expression() const;

 private:
  void init(const std::wstring& string);

  /// convert wchar to index
  Type wchar_to_index(const wchar_t);

  /// convet wchar with prime, for example a', to index
  Type wchar_with_prime_to_index(const wchar_t);

 private:
  Type index_;
  Spin spin_;

  /// the name that user passed in from the constructor
  std::wstring name_;
};

template <typename String, typename>
OrbitalIndex::OrbitalIndex(String &&symbol) {

  name_ = utility::to_wstring(symbol);

  if (name_.find_first_of(L'_') == std::wstring::npos) {
    init(name_);
    spin_ = Spin::None;
  } else {
    auto left = name_.find_first_of(L'_');

    TA_ASSERT(left != std::wstring::npos);
    TA_ASSERT(name_.size() - left == 2);
    std::wstring sub_letter = name_.substr(0, left);

    init(sub_letter);

    wchar_t spin = name_[left + 1];

    if (spin == L'α') {
      spin_ = Spin::Alpha;
    } else if (spin == L'β') {
      spin_ = Spin::Beta;
    } else {
      throw std::runtime_error("Wrong Spin Label");
    }
  }
}

}

#endif  // MPQC_ORBITAL_INDEX_H
