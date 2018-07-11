
#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_BASIS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_BASIS_H_

#include <iosfwd>
#include <memory>
#include <vector>

#include <libint2/basis.h>
#include <madness/world/array_addons.h>
#include <tiledarray.h>
#include "mpqc/util/external/boost/small_vector.h"

#include "mpqc/util/keyval/keyval.h"

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/util/misc/observer.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

using Shell = libint2::Shell;
using ShellVec = std::vector<Shell>;

/// @ingroup ChemistryESLCAOBasis
/// @{

/// Basis is a clustered sequence of libint2::Shell objects.
/// The sequence is represented as a vector of vectors of shells.
class Basis : virtual public DescribedClass {
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

    /// The KeyVal ctor queries the following keywords:
    /// | Keyword | Type | Default| Description |
    /// |---------|------|--------|-------------|
    ///  |\c name |string|none    |The basis set name|
    ///
    Factory(const KeyVal& kv);

    /// This ctor creates a Factory that will use default_name for all atoms
    Factory(std::string const & default_name);

    /*! \brief Constructs a vector of ShellVecs
     *
     * Each ShellVec represents the shells for a cluster.
     */
    std::vector<ShellVec> get_cluster_shells(Molecule const &) const;

    /*! \brief returns a single vector of all shells in the molecule */
    ShellVec get_flat_shells(Molecule const &) const;

    /// @return the basis set name
    const std::string& name() const { return basis_set_name_; }

   private:
    std::string basis_set_name_;
  };  // Basis::Factory

  /// creates an empty Basis
  Basis();
  ~Basis();
  Basis(Basis const &);
  Basis(Basis &&);
  Basis &operator=(Basis const &);
  Basis &operator=(Basis &&);

  /// construct a Basis from a clustered sequence of shells
  Basis(std::vector<ShellVec> shells) : shells_(std::move(shells)) {}

  /// @return a reference to the vector of shell clusters
  std::vector<ShellVec> const &cluster_shells() const;
  /// @return the vector of all shells
  std::vector<Shell> flattened_shells() const;

  TiledArray::TiledRange1 create_trange1() const;

  /// @return the maximum number of primitives in any shell in this Basis
  int64_t max_nprim() const;
  /// @return the highest angular momentum of any shell in this Basis
  int64_t max_am() const;
  /// @return the total number of functions in this Basis
  int64_t nfunctions() const;
  /// @return the total number of shells in this Basis
  int64_t nshells() const;
  /// @return the total number of shell clusters in this Basis
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

    for (auto i = 0ul; i < nvecs; ++i) {
      ShellVec tmp;
      ar &tmp;
      shells_.emplace_back(std::move(tmp));
    }
  }

 protected:
  std::vector<ShellVec> shells_;

  friend Basis merge(const Basis &basis1, const Basis &basis2);
};

/// merges two Basis objects by concatenating their shell cluster sequences.
/// @return a Basis object in which shells of \c basis1 and \c basis2
Basis merge(const Basis &basis1, const Basis &basis2);

std::ostream& operator<<(std::ostream &os, Basis::Factory const &f);

std::ostream& operator<<(std::ostream &, Basis const &);

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
 * construct a map that maps column of basis to column of sub_basis, sub_basis has to be a subset of basis
 * @warning the value in index starts with 1, value 0 in index indicates this column is missing in sub_basis
 * This approach uses N^2 algorithm
 * //TODO need unit test for this
 * @param basis Basis object
 * @param sub_basis Basis object, subset of basis
 * @return a vector of column id of sub_basis
 */
Eigen::RowVectorXi sub_basis_map(const Basis& basis, const Basis& sub_basis);

/**
 * construct Basis from a factory and a Molecule on process 0 and broadcast to the entire world
 * @param world the madness::World
 * @param factory the Basis::Factory object
 * @param mol the Molecule object
 */
Basis parallel_make_basis(madness::World &world, const Basis::Factory& factory,
                          const mpqc::Molecule &mol);

/// AtomicBasis is a Basis constructed from a collection of atoms, i.e.
/// a Molecule( or a UnitCell). Updating
/// the Molecule object updates this basis as well. Clustering of AtomicBasis
/// does not have to correspond to the cluster structure of Molecule (e.g.
/// as a result of reblocking).
/// @note Shells are organized in the order of increasing atom index.
class AtomicBasis : public Basis, public utility::Observer {
 public:
  using Shell = ::mpqc::lcao::gaussian::Shell;

  /**
   * \brief the KeyVal constructor for MolecularBasis
   *
   * The KeyVal constructor queries all keywords of Basis and Basis::Factory, as well
   * as the following additional keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * |atoms| Molecule or UnitCell | none | the Molecule object |
   * |reblock| int | 0 (i.e., no reblocking) | size to reblock basis |
   */
  AtomicBasis(const KeyVal& kv);

  std::shared_ptr<const Factory> factory() const;
  std::shared_ptr<const Molecule> molecule() const;

  /// @return the string identifier for this basis, i.e. factory()->name()
  const std::string& name() const { return factory()->name(); }

 private:
  std::shared_ptr<Factory> factory_;
  std::shared_ptr<Molecule> molecule_;
  std::vector<std::vector<size_t>> shell_to_atom_;  // maps shells_ to the atoms in molecule_

  void compute_shell_to_atom();
  // updates shells with the atomic coordinates
  void rebuild_shells();
};

std::ostream &operator<<(std::ostream &, AtomicBasis const &);

/// @}

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
#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_BASIS_H_
