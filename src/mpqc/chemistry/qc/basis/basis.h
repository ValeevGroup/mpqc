
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_H_

#include <iosfwd>
#include <memory>
#include <vector>

#include <libint2/basis.h>
#include <madness/world/array_addons.h>
#include <tiledarray.h>

#include "mpqc/util/keyval/keyval.h"

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/basis/basis_fwd.h"
#include "mpqc/util/misc/observer.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

using Shell = libint2::Shell;
using ShellVec = std::vector<Shell>;

/*
 * \defgroup Basis Basis
 *
 * \brief The Basis module contains information about basis set
 *
 *
 */

/// Basis is a clustered sequence of libint2::Shell objects.
/// The sequence is represented as a vector of vectors of shells.
class Basis : virtual public DescribedClass,
              public utility::Observer {
 public:
  using Shell = ::mpqc::lcao::gaussian::Shell;

  /// Factory is a ctor helper
  class Factory {
   public:
    Factory() = delete;  // Can't init a basis without name.
    Factory(Factory const &b) = default;
    Factory(Factory &&b) = default;
    Factory &operator=(Factory const &b) = default;
    Factory &operator=(Factory &&b) = default;

    /// This ctor creates a Factory that will use default_name for all atoms
    Factory(std::string const & default_name);

    /*! \brief Constructs a vector of ShellVecs
     *
     * Each ShellVec represents the shells for a cluster.
     */
    std::vector<ShellVec> get_cluster_shells(Molecule const &) const;

    /*! \brief returns a single vector of all shells in the molecule */
    ShellVec get_flat_shells(Molecule const &) const;

   private:
    std::string basis_set_name_;
  };

  /// created an empty Basis
  Basis();
  ~Basis();
  Basis(Basis const &);
  Basis(Basis &&);
  Basis &operator=(Basis const &);
  Basis &operator=(Basis &&);

  /// constructs a Basis object from a vector of shell clusters
  explicit Basis(std::vector<ShellVec> cs);

  /**
   * \brief the KeyVal constructor for Basis
   *
   * The KeyVal constructor uses the following keywords:
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * |molecule|Molecule|none|keyval to construct molecule|
   * |name|string|none|basis set name|
   * |reblock| int | 0 | size to reblock basis |
   *
   */
  Basis(const KeyVal &kv);

  /// @return a reference to the vector of shell clusters
  std::vector<ShellVec> const &cluster_shells() const;
  /// @return the vector of all shells
  std::vector<Shell> flattened_shells() const;

  TiledArray::TiledRange1 create_trange1() const;

  /// @return the maximum number of primitives in any shell in this Basis
  int64_t max_nprim() const;
  /// @return the highest angular momentum of any shell in this Basis
  int64_t max_am() const;
  /// @return the maximum number of functions (i.e., size) of any shell in this Basis
  int64_t nfunctions() const;
  /// @return the number of shells in this Basis
  int64_t nshells() const;
  /// @return the number of shell clusters in this Basis
  int64_t nclusters() const { return shells_.size(); };

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

  friend Basis merge(const Basis &basis1, const Basis &basis2);
};

/// merges two Basis objects by concatenating their shell cluster sequences.
/// @return a Basis object in which shells of \c basis1 and \c basis2
Basis merge(const Basis &basis1, const Basis &basis2);

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

/**
 * construct basis on MPI process 0 and broadcast to all processes
 * @param world madness::World
 * @param basis_factory Basis::Factory object
 * @param mol Molecule object
 * @return Basis object
 */
Basis parallel_construct_basis(madness::World &world, const Basis::Factory &basis_factory,
                               const mpqc::Molecule &mol);
/**
 * construct a map that maps column of basis to column of sub_basis, sub_basis has to be a subset of basis
 * @warning the value in index starts with 1, value 0 in index indicates this column is missing in sub_basis
 * This approach uses N^2 algorithm
 * //TODO need unit test for this
 * @param basis Basis object
 * @param sub_basis Basis object, subset of basis
 * @return a vector of column id of sub_basis
 */
Eigen::RowVectorXi sub_basis_map(const Basis& basis, const Basis& sub_basis);

}  // namespace gaussian
}  // namespace lcao
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

}  // namespace madness
}  // namespace archive
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_H_
