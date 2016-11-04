
#ifndef MPQC_BASIS_BASIS_H
#define MPQC_BASIS_BASIS_H

#include <iosfwd>
#include <memory>
#include <vector>

#include <libint2/shell.h>
#include <madness/world/array_addons.h>
#include <tiledarray.h>

#include "mpqc/util/keyval/keyval.h"

#include "mpqc/chemistry/molecule/molecule_fwd.h"
#include "mpqc/chemistry/qc/basis/basis_fwd.h"
#include "mpqc/chemistry/qc/basis/basis_set.h"

namespace mpqc {
namespace basis {

using Shell = libint2::Shell;
using ShellVec = std::vector<Shell>;

/*
 * \defgroup Basis Basis
 *
 * \brief The Basis module contains information about basis set
 *
 *
 */

class Basis : public DescribedClass {
 public:
  using Shell = libint2::Shell;

  Basis();
  ~Basis();
  Basis(Basis const &);
  Basis(Basis &&);
  Basis &operator=(Basis const &);
  Basis &operator=(Basis &&);

  ///
  Basis(std::vector<ShellVec> cs);

  /**
   * \brief KeyVal constructor for Basis
   *
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * |molecule|Molecule|none|keyval to construct molecule|
   * |name|string|none|basis set name|
   *
   */
  Basis(const KeyVal &kv);

  /// join another basis together
  Basis join(const Basis &basis);

  std::vector<ShellVec> const &cluster_shells() const;

  TiledArray::TiledRange1 create_trange1() const;

  int64_t max_nprim() const;
  int64_t max_am() const;
  int64_t nfunctions() const;
  int64_t nshells() const;
  int64_t nclusters() const { return shells_.size(); };
  std::vector<Shell> flattened_shells() const;

  template <typename Archive>
  typename std::enable_if<
      madness::archive::is_output_archive<Archive>::value>::type
  serialize(const Archive &ar) {
    std::size_t nvecs = shells_.size();
    ar &nvecs;
    for (auto const &v : shells_) {
      ar &v;
    }
  }

  template <typename Archive>
  typename std::enable_if<
      madness::archive::is_input_archive<Archive>::value>::type
  serialize(const Archive &ar) {
    std::size_t nvecs = 0;
    ar &nvecs;

    for (auto i = 0; i < nvecs; ++i) {
      ShellVec tmp;
      ar &tmp;
      shells_.emplace_back(std::move(tmp));
    }
  }

 private:
  std::vector<ShellVec> shells_;
};

std::ostream &operator<<(std::ostream &, Basis const &);

/*! \brief reblock allows for reblocking a basis
 *
 * \warning If reblocking a basis with the intent to use it with tensors
 * computed with the old basis you must be careful not to reorder the shells.
 *
 * \param op should be a function that takes a std::vector<Shell> and returns
 * a std::vector<std::vector<Shell>> for use in initializing a Basis.
 */
template <typename Op, typename... Args>
Basis reblock(Basis const &basis, Op op, Args... args) {
  return Basis(op(basis.flattened_shells(), args...));
}

Basis parallel_construct_basis(madness::World &world, const BasisSet &basis_set,
                               const mpqc::Molecule &mol);

}  // namespace basis
}  // namespace mpqc

namespace madness {
namespace archive {

// serialize libint2::Shell object

template <typename Archive>
struct ArchiveSerializeImpl<Archive, libint2::Shell::Contraction> {
  static inline void serialize(const Archive &ar,
                               libint2::Shell::Contraction &c) {
    ar &c.l;
    ar &c.pure;
    ar &c.coeff;
  }
};

template <typename Archive>
struct ArchiveSerializeImpl<Archive, libint2::Shell> {
  static inline void serialize(const Archive &ar, libint2::Shell &s) {
    ar &s.alpha;
    ar &s.contr;
    ar &s.O;
    ar &s.max_ln_coeff;
  };
};

}  // end of namespace madness
}  // end of namespace archive
#endif /* end of include guard: MPQC_BASIS_BASIS_H */
