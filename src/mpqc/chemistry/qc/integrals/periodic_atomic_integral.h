#ifndef MPQC_PERIODIC_ATOMIC_INTEGRAL_H
#define MPQC_PERIODIC_ATOMIC_INTEGRAL_H

#include <iosfwd>
#include <vector>

#include "atomic_integral_base.h"
#include "atomic_integral_base.cpp"
#include "atomic_integral.h"
#include <mpqc/chemistry/qc/integrals/integrals.h>
#include <mpqc/util/keyval/keyval.hpp>
// Eigen matrix algebra library
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::Vector3i Vec3I;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor> Matrixc;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Vectorc;

// constant
const std::complex<double> I(0.0, 1.0);
const auto angstrom_to_bohr = 1 / 0.52917721092;  // 2010 CODATA value

using namespace mpqc::integrals::detail;

namespace mpqc {
namespace integrals {

template <typename Tile, typename Policy>
class PeriodicAtomicIntegral : public AtomicIntegralBase {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  //    using Op = std::function<Tile(TA::TensorD&&)>;
  using Op = std::function<Tile(TA::TensorZ &&)>;

  PeriodicAtomicIntegral() = default;
  PeriodicAtomicIntegral(PeriodicAtomicIntegral &&) = default;
  PeriodicAtomicIntegral &operator=(PeriodicAtomicIntegral &&) = default;

  PeriodicAtomicIntegral(const KeyVal &kv)
      : AtomicIntegralBase(kv),
        ao_formula_registry_(),
        orbital_space_registry_() {
    std::string molecule_type = kv.value<std::string>("molecule:type");

    if (molecule_type != "PeriodicSystem") {
      throw std::invalid_argument("Moleule Type Has To be PeriodicSystem!!");
    }
    auto ps = std::dynamic_pointer_cast<molecule::PeriodicSystem>(mol_);
    dcell_ = ps->dcell();
    R_max_ = ps->R_max();
    RD_max_ = ps->RD_max();
    RJ_max_ = ps->RJ_max();
    nk_ = ps->nk();
    op_ = mpqc::ta_routines::TensorZPassThrough();
  }

  ~PeriodicAtomicIntegral() noexcept = default;

  /// wrapper to compute function
  //  std::vector<TArray> compute(const std::wstring &);

  /**
   *  compute integral by Formula
   *  this function will look into registry first
   *  if Formula computed, it will return it from registry
   *  if not, it will compute it
   */
  //  std::vector<TArray> compute(const Formula &formula);

  /// wrapper to compute function
  TArray compute(const std::wstring &);

  /**
   *  compute integral by Formula
   *  this function will look into registry first
   *  if Formula computed, it will return it from registry
   *  if not, it will compute it
   */
  TArray compute(const Formula &formula);

  /// transform a matrix from real to reciprocal space
  TArray transform_real2recip(TArray &matrix);

  /// compute density: D = C(occ).C(occ)t
  TArray compute_density(TArray &fock_real, TArray &fock_recip,
                         TArray &overlap,
                         int64_t ndocc);


 private:
  FormulaRegistry<TArray> ao_formula_registry_;
  std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
  Op op_;

  // Density
  std::vector<TArray> D_;
  TArray density_;

  // Constants

  Vec3I R_max_ = {0, 0,
                  0};  // range of expansion of Bloch Gaussians in AO Gaussians
  Vec3I RJ_max_ = {0, 0, 0};       // range of Coulomb operation
  Vec3I RD_max_ = {0, 0, 0};       // range of density representation
  Vec3I nk_ = {1, 1, 1};           // # of k points in each direction
  Vec3D dcell_ = {0.0, 0.0, 0.0};  // direct unit cell params (in a.u.)

  int64_t idx_lattice(int x, int y, int z, Vec3I vec);
  Vec3D R_vector(int64_t idx_lattice, Vec3I vec);
  int64_t idx_k(int x, int y, int z, Vec3I nk);
  Vec3D k_vector(int64_t idx_k);

  template <typename U = Policy>
  TA::DistArray<
      Tile, typename std::enable_if<std::is_same<U, TA::SparsePolicy>::value,
                                    TA::SparsePolicy>::type>
  compute_integrals(
      madness::World &world, ShrPool<libint2::Engine> &engine,
      Bvector const &bases,
      std::shared_ptr<Screener> p_screen =
          std::make_shared<integrals::Screener>(integrals::Screener{})) {
    auto result = sparse_complex_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

  void parse_one_body_periodic(
      const Formula &formula,
      std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool, Bvector &bases,
      molecule::Molecule &shifted_mol);

  void parse_two_body_periodic(const Formula &formula,
      std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool, Bvector &bases,
      Vec3D shift_coul, bool if_coulomb);

  std::shared_ptr<basis::Basis> shift_basis_origin(basis::Basis &basis,
                                                   Vec3D shift);

  std::shared_ptr<basis::Basis> shift_basis_origin(basis::Basis &basis,
                                                   Vec3D shift_base,
                                                   Vec3I nshift,
                                                   bool is_real_space);

  std::shared_ptr<molecule::Molecule> shift_mol_origin(molecule::Molecule &mol,
                                                       Vec3D shift);

  libint2::any to_libint2_operator_params(Operator::Type mpqc_oper,
                                          const AtomicIntegralBase &base,
                                          molecule::Molecule &mol);

  template <typename E>
  TA::DistArray<Tile, TA::SparsePolicy> sparse_complex_integrals(
      mad::World &world, ShrPool<E> shr_pool, Bvector const &bases,
      std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
      std::function<Tile(TA::TensorZ &&)> op =
          mpqc::ta_routines::TensorZPassThrough());

  void sort_eigen(Vectorc &eigVal, Matrixc &eigVec);
};

template <typename Tile, typename Policy>
typename PeriodicAtomicIntegral<Tile, Policy>::TArray PeriodicAtomicIntegral<
    Tile, Policy>::compute(const std::wstring &formula_string) {
  auto formula = Formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename PeriodicAtomicIntegral<Tile, Policy>::TArray
PeriodicAtomicIntegral<Tile, Policy>::compute(const Formula &formula) {

  TArray result;
  Bvector bs_array;
  std::shared_ptr<EnginePool<libint2::Engine>> engine_pool;
  double size = 0.0;

  if (formula.rank() == 2) {
    if (formula.oper().type() == Operator::Type::Kinetic
            || formula.oper().type() == Operator::Type::Overlap) {
      parse_one_body_periodic(formula, engine_pool, bs_array, *mol_);
      result = compute_integrals(this->world_, engine_pool, bs_array);
    }
    else if (formula.oper().type() == Operator::Type::Nuclear) {

      auto RJ_size =
          1 + idx_lattice(RJ_max_(0), RJ_max_(1), RJ_max_(2), RJ_max_);
      for (auto RJ = 0; RJ < RJ_size; ++RJ) {
          auto shift_mol = R_vector(RJ, RJ_max_);
          auto shifted_mol = shift_mol_origin(*mol_, shift_mol);
          parse_one_body_periodic(formula, engine_pool, bs_array,
                                  *shifted_mol);
          if (RJ == 0)
              result = compute_integrals(this->world_, engine_pool, bs_array);
          else
              result("mu, nu") +=
                      compute_integrals(this->world_, engine_pool, bs_array)
                      ("mu, nu");
      }
    }
    else
        throw std::runtime_error("Rank-2 operator type not supported");

    size = utility::array_size(result);
    utility::print_par(world_,
                       "\nComputed One Body Integral for Periodic System: ",
                       utility::to_string(formula.string()));
    utility::print_par(world_, " Size: ", size, " GB\n");

  }
  else if (formula.rank() == 4) {
      if (formula.oper().type() == Operator::Type::J) {
          auto j_formula = formula;
          j_formula.set_operator_type(Operator::Type::Coulomb);
          auto RJ_size = 1 + idx_lattice(RJ_max_(0), RJ_max_(1), RJ_max_(2), RJ_max_);
          for (auto RJ = 0; RJ < RJ_size; ++RJ) {
              auto vec_RJ = R_vector(RJ, RJ_max_);
              parse_two_body_periodic(j_formula, engine_pool, bs_array, vec_RJ, true);
              auto J = compute_integrals(this->world_, engine_pool, bs_array);
              if (RJ == 0)
                  result("mu, nu") = J("mu, nu, lambda, rho") * density_("lambda, rho");
              else
                  result("mu, nu") += J("mu, nu, lambda, rho") * density_("lambda, rho");
          }
      }
      else if (formula.oper().type() == Operator::Type::K) {
          auto k_formula = formula;
          k_formula.set_operator_type(Operator::Type::Coulomb);
          auto RJ_size = 1 + idx_lattice(RJ_max_(0), RJ_max_(1), RJ_max_(2), RJ_max_);
          for (auto RJ = 0; RJ < RJ_size; ++RJ) {
              auto vec_RJ = R_vector(RJ, RJ_max_);
              parse_two_body_periodic(k_formula, engine_pool, bs_array, vec_RJ, false);
              auto K = compute_integrals(this->world_, engine_pool, bs_array);
              if (RJ == 0)
                  result("mu, nu") = K("mu, lambda, nu, rho") * density_("lambda, rho");
              else
                  result("mu, nu") += K("mu, lambda, nu, rho") * density_("lambda, rho");
          }
      }
      else
          throw std::runtime_error("Rank-4 operator type not supported");

  }
  else
      throw std::runtime_error("Operator rank not supported");

  return result;
}


template <typename Tile, typename Policy>
void PeriodicAtomicIntegral<Tile, Policy>::parse_one_body_periodic(
        const Formula &formula,
        std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool, Bvector &bases,
        molecule::Molecule &shifted_mol) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 1);
  TA_ASSERT(ket_indices.size() == 1);

  auto bra_index = bra_indices[0];
  auto ket_index = ket_indices[0];

  TA_ASSERT(bra_index.is_ao());
  TA_ASSERT(ket_index.is_ao());

  auto bra_basis = this->index_to_basis(bra_index);
  auto ket_basis = this->index_to_basis(ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  // Form a compound ket basis by shifting origins from -Rmax to Rmax
  Vec3D zero_shift_base(0.0, 0.0, 0.0);
  ket_basis = shift_basis_origin(*ket_basis, zero_shift_base, R_max_, true);

  bases = Bvector{{*bra_basis, *ket_basis}};

  auto oper_type = formula.oper().type();
  engine_pool = integrals::make_engine_pool(
      to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis, *ket_basis), libint2::BraKet::x_x,
      to_libint2_operator_params(oper_type, *this, shifted_mol));
}

template <typename Tile, typename Policy>
void PeriodicAtomicIntegral<Tile, Policy>::parse_two_body_periodic(const Formula &formula,
    std::shared_ptr<EnginePool<libint2::Engine>> &engine_pool, Bvector &bases,
    Vec3D shift_coul, bool if_coulomb) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 2);
  TA_ASSERT(ket_indices.size() == 2);

  auto bra_index0 = bra_indices[0];
  auto bra_index1 = bra_indices[0];
  auto ket_index0 = ket_indices[0];
  auto ket_index1 = ket_indices[0];

  TA_ASSERT(bra_index0.is_ao());
  TA_ASSERT(bra_index1.is_ao());
  TA_ASSERT(ket_index0.is_ao());
  TA_ASSERT(ket_index1.is_ao());

  auto bra_basis0 = this->index_to_basis(bra_index0);
  auto bra_basis1 = this->index_to_basis(bra_index1);
  auto ket_basis0 = this->index_to_basis(ket_index0);
  auto ket_basis1 = this->index_to_basis(ket_index1);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(bra_basis1 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  // Form a compound index basis
  Vec3D zero_shift_base(0.0, 0.0, 0.0);
  if (if_coulomb) {
      bra_basis1 = shift_basis_origin(*bra_basis1, zero_shift_base, R_max_, true);
      ket_basis0 = shift_basis_origin(*ket_basis0, shift_coul);
  }
  else {
      bra_basis1 = shift_basis_origin(*bra_basis1, shift_coul);
      ket_basis0 = shift_basis_origin(*ket_basis0, zero_shift_base, R_max_, true);
  }
  ket_basis1 = shift_basis_origin(*ket_basis1, shift_coul, RD_max_, true);

  if (formula.notation() == Formula::Notation::Chemical)
    bases = Bvector{{*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1}};
  else
      throw "Physical notation not supported!";
//    bases = Bvector{{*bra_basis0, *ket_basis0, *bra_basis1, *ket_basis1}};

  auto oper_type = formula.oper().type();
  engine_pool = integrals::make_engine_pool(
      to_libint2_operator(oper_type),
      utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
      libint2::BraKet::xx_xx,
      to_libint2_operator_params(oper_type, *this, *mol_));
}

template <typename Tile, typename Policy>
int64_t PeriodicAtomicIntegral<Tile, Policy>::idx_lattice(int x, int y, int z,
                                                          Vec3I vec) {
  if (vec(0) >= 0 && vec(1) >= 0 && vec(2) >= 0 && abs(x) <= vec(0) &&
      abs(y) <= vec(1) && abs(z) <= vec(2)) {
    int64_t idx = (x + vec(0)) * (2 * vec(0) + 1) * (2 * vec(1) + 1) +
                  (y + vec(1)) * (2 * vec(1) + 1) + (z + vec(2));
    return idx;
  } else {
    throw "invalid lattice sum index/boundaries";
  }
}

template <typename Tile, typename Policy>
int64_t PeriodicAtomicIntegral<Tile, Policy>::idx_k(int x, int y, int z,
                                                    Vec3I nk) {
  if (nk(0) >= 1 && nk(1) >= 1 && nk(2) >= 1 && x >= 0 && y >= 0 && z >= 0 &&
      x < nk(0) && y < nk(1) && z < nk(2)) {
    int64_t idx = x * nk(0) * nk(1) + y * nk(1) + z;
    return idx;
  } else {
    throw "invalid k-space index/boundaries";
  }
}

template <typename Tile, typename Policy>
Vec3D PeriodicAtomicIntegral<Tile, Policy>::k_vector(int64_t idx_k) {
  Vec3D result;
  auto x = idx_k / nk_(2) / nk_(1);
  auto y = (idx_k / nk_(2)) % nk_(1);
  auto z = idx_k % nk_(2);
  result(0) =
      (dcell_(0) == 0.0) ? 0.0 : (-1.0 + (2.0 * (x + 1) - 1.0) / nk_(0)) *
                                     (M_PI / dcell_(0));
  result(1) =
      (dcell_(1) == 0.0) ? 0.0 : (-1.0 + (2.0 * (y + 1) - 1.0) / nk_(1)) *
                                     (M_PI / dcell_(1));
  result(2) =
      (dcell_(2) == 0.0) ? 0.0 : (-1.0 + (2.0 * (z + 1) - 1.0) / nk_(2)) *
                                     (M_PI / dcell_(2));
  return result;
}

template <typename Tile, typename Policy>
Vec3D PeriodicAtomicIntegral<Tile, Policy>::R_vector(int64_t idx_lattice,
                                                     Vec3I vec) {
  auto z = idx_lattice % (2 * vec(2) + 1);
  auto y = (idx_lattice / (2 * vec(2) + 1)) % (2 * vec(1) + 1);
  auto x = idx_lattice / (2 * vec(2) + 1) / (2 * vec(1) + 1);
  Vec3D result((x - vec(0)) * dcell_(0), (y - vec(1)) * dcell_(1),
               (z - vec(2)) * dcell_(2));
  return result;
}

template <typename Tile, typename Policy>
std::shared_ptr<basis::Basis>
PeriodicAtomicIntegral<Tile, Policy>::shift_basis_origin(basis::Basis &basis,
                                                         Vec3D shift) {
  std::vector<ShellVec> vec_of_shells;
  for (auto shell_vec : basis.cluster_shells()) {
    ShellVec shells;
    for (auto shell : shell_vec) {
      std::array<double, 3> new_origin = {
          shell.O[0] + shift(0), shell.O[1] + shift(1), shell.O[2] + shift(2)};
      shell.move(new_origin);
      shells.push_back(shell);
    }
    vec_of_shells.push_back(shells);
  }
  basis::Basis result(vec_of_shells);
  auto result_ptr = std::make_shared<basis::Basis>(result);
  return result_ptr;
}

// Form a compound basis with double index (R, u):
// R is cell index : [-R_max, R_max]
// u is conventional AO index within a cell
template <typename Tile, typename Policy>
std::shared_ptr<basis::Basis>
PeriodicAtomicIntegral<Tile, Policy>::shift_basis_origin(basis::Basis &basis,
                                                         Vec3D shift_base,
                                                         Vec3I nshift,
                                                         bool is_real_space) {
  std::vector<ShellVec> vec_of_shells;

  int64_t shift_size =
      is_real_space ? (1 + idx_lattice(nshift(0), nshift(1), nshift(2), nshift))
                    : (1 + idx_k(nshift(0), nshift(1), nshift(2), nshift));

  for (auto idx_shift = 0; idx_shift < shift_size; ++idx_shift) {
    Vec3D shift = is_real_space ? (R_vector(idx_shift, nshift) + shift_base)
                                : (k_vector(idx_shift) + shift_base);

    for (auto shell_vec : basis.cluster_shells()) {
      ShellVec shells;
      for (auto shell : shell_vec) {
        std::array<double, 3> new_origin = {shell.O[0] + shift(0),
                                            shell.O[1] + shift(1),
                                            shell.O[2] + shift(2)};
        shell.move(new_origin);
        shells.push_back(shell);
      }
      vec_of_shells.push_back(shells);
    }
  }

  basis::Basis result(vec_of_shells);
  auto result_ptr = std::make_shared<basis::Basis>(result);
  return result_ptr;
}

template <typename Tile, typename Policy>
std::shared_ptr<molecule::Molecule>
PeriodicAtomicIntegral<Tile, Policy>::shift_mol_origin(molecule::Molecule &mol,
                                                       Vec3D shift) {
  std::vector<molecule::AtomBasedClusterable> vec_of_clusters;
  for (auto &cluster : mol) {
    molecule::AtomBasedCluster shifted_cluster;
    for (auto &atom : molecule::collapse_to_atoms(cluster)) {
      molecule::Atom shifted_atom(atom.center() + shift, atom.mass(),
                                  atom.charge());
      shifted_cluster.add_clusterable(shifted_atom);
    }
    shifted_cluster.update_cluster();
    vec_of_clusters.emplace_back(shifted_cluster);
  }

  molecule::Molecule result(vec_of_clusters);

  auto result_ptr = std::make_shared<molecule::Molecule>(result);
  return result_ptr;
}

template <typename Tile, typename Policy>
libint2::any PeriodicAtomicIntegral<Tile, Policy>::to_libint2_operator_params(
    Operator::Type mpqc_oper, const AtomicIntegralBase &base,
    molecule::Molecule &mol) {
  TA_USER_ASSERT((Operator::Type::__first_1body_operator <= mpqc_oper &&
                  mpqc_oper <= Operator::Type::__last_1body_operator) ||
                     (Operator::Type::__first_2body_operator <= mpqc_oper &&
                      mpqc_oper <= Operator::Type::__last_2body_operator),
                 "invalid Operator::Type");

  libint2::any result;
  switch (mpqc_oper) {
    case Operator::Type::Nuclear: {
      result = make_q(mol);
    } break;
    case Operator::Type::cGTG:
    case Operator::Type::cGTGCoulomb:
    case Operator::Type::DelcGTG2: {
      result = base.gtg_params();
    } break;
    case Operator::Type::cGTG2: {
      const auto &cgtg_params = base.gtg_params();
      const auto ng = cgtg_params.size();
      std::decay<decltype(cgtg_params)>::type cgtg2_params;
      cgtg2_params.reserve(ng * (ng + 1) / 2);
      for (auto b = 0; b < ng; ++b) {
        for (auto k = 0; k <= b; ++k) {
          const auto gexp = cgtg_params[b].first + cgtg_params[k].first;
          const auto gcoeff = cgtg_params[b].second * cgtg_params[k].second *
                              (b == k ? 1 : 2);  // if a != b include ab and ba
          cgtg2_params.push_back(std::make_pair(gexp, gcoeff));
        }
      }
      result = cgtg2_params;
    } break;
    default:;  // nothing to do
  }
  return result;
}

/*! \brief Construct sparse complex integral tensors in parallel.
 *
 * \param shr_pool should be a std::shared_ptr to an IntegralEnginePool
 * \param bases should be a std::array of Basis, which will be copied.
 * \param op needs to be a function or functor that takes a TA::TensorZ && and
 * returns any valid tile type. Op is copied so it can be moved.
 * ```
 * auto t = [](TA::TensorZ &&ten){return std::move(ten);};
 * ```
 *
 * \param screen should be a std::shared_ptr to a Screener.
 */
template <typename Tile, typename Policy>
template <typename E>
TA::DistArray<Tile, TA::SparsePolicy>
PeriodicAtomicIntegral<Tile, Policy>::sparse_complex_integrals(
    mad::World &world, ShrPool<E> shr_pool, Bvector const &bases,
    std::shared_ptr<Screener> screen, std::function<Tile(TA::TensorZ &&)> op) {
  // Build the Trange and Shape Tensor
  auto trange = detail::create_trange(bases);
  const auto tvolume = trange.tiles_range().volume();
  std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
  TA::TensorF tile_norms(trange.tiles_range(), 0.0);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<Bvector>(bases);

  // Make a pointer to an Integral builder.  Doing this because we want to use
  // it in Tasks.
  auto builder_ptr = std::make_shared<IntegralBuilder<Tile, E>>(
      std::move(shr_pool), std::move(shr_bases), std::move(screen),
      std::move(op));

  auto task_f = [=](int64_t ord, detail::IdxVec idx, TA::Range rng,
                    TA::TensorF *tile_norms_ptr, Tile *out_tile) {

    // This is why builder was made into a shared_ptr.
    auto &builder = *builder_ptr;
    auto ta_tile = builder.integrals(idx, std::move(rng));

    const auto tile_volume = ta_tile.range().volume();
    const auto tile_norm = ta_tile.norm();

    // Keep tile if it was significant.
    bool save_norm =
        tile_norm >= tile_volume * TA::SparseShape<float>::threshold();
    if (save_norm) {
      *out_tile = builder.op(std::move(ta_tile));

      auto &norms = *tile_norms_ptr;
      norms[ord] = tile_norm;
    }
  };

  auto pmap = SpPolicy::default_pmap(world, tvolume);
  for (auto const ord : *pmap) {
    tiles[ord].first = ord;
    detail::IdxVec idx = trange.tiles_range().idx(ord);
    world.taskq.add(task_f, ord, idx, trange.make_tile_range(ord), &tile_norms,
                    &tiles[ord].second);
  }
  world.gop.fence();

  TA::SparseShape<float> shape(world, tile_norms, trange);
  TA::DistArray<Tile, SpPolicy> out(world, trange, shape, pmap);

  detail::set_array(tiles, out);
  out.truncate();

  return out;
}

template <typename Tile, typename Policy>
typename PeriodicAtomicIntegral<Tile, Policy>::TArray
PeriodicAtomicIntegral<Tile, Policy>::transform_real2recip(
    TArray &matrix) {
  TArray result;
  auto k_size = 1 + idx_k(nk_(0) - 1, nk_(1) - 1, nk_(2) - 1, nk_);
  auto R_size = 1 + idx_lattice(R_max_(0), R_max_(1), R_max_(2), R_max_);
  auto tr0 = matrix.trange().data()[0];

  // Make tiled range for the compound index (k,v)
  auto blocking = std::vector<int64_t> {0};
  for (auto k = 0; k < k_size; ++k) {
      for (auto u = 0; u < tr0.tile_extent(); ++u) {
          auto next = blocking.back() + tr0.tile(u).second - tr0.tile(u).first;
          blocking.emplace_back(next);
      }
  }
  TA::TiledRange1 tr1(blocking.begin(), blocking.end());

  // Perform real->reciprocal transformation with Eigen
  // TODO: perform it with TA (take arg tile from "matrix",
  // transform it, add it to result tile in "result", construct pmap&shape,
  // use MADNESSworld ...)

  auto matrix_eig = array_ops::array_to_eigen(matrix);
  Eig::MatrixXcd result_eig(tr0.extent(), tr1.extent());
  result_eig.setZero();

  auto threshold = std::numeric_limits<double>::epsilon();
  for (auto R = 0; R < R_size; ++R) {
      auto bmat = matrix_eig.block(0, R*tr0.extent(), tr0.extent(), tr0.extent());
      if (bmat.norm() < bmat.size() * threshold)
          continue;
      else {
          auto vec_R = R_vector(R, R_max_);
          for (auto k = 0; k < k_size; ++k) {
              auto vec_k = k_vector(k);
              auto exponent = std::exp(I * vec_k.dot(vec_R));
              result_eig.block(0, k*tr0.extent(), tr0.extent(), tr0.extent()) +=
                      bmat * exponent;
          }
      }
  }

  result = array_ops::eigen_to_array<Tile>(world_, result_eig, tr0, tr1);

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAtomicIntegral<Tile, Policy>::TArray
PeriodicAtomicIntegral<Tile, Policy>::compute_density(
    TArray &fock_real, TArray &fock_recip,
    TArray &overlap,
    int64_t ndocc) {

    TArray result;

    auto k_size = 1 + idx_k(nk_(0) - 1, nk_(1) - 1, nk_(2) - 1, nk_);

    std::vector<Vectorc> eps(k_size);
    std::vector<Matrixc> C(k_size);

    auto tr0 = fock_real.trange().data()[0];
    //TODO: write a function to form tr1. Now R_max_ must equal RD_max_
    auto tr1 = fock_real.trange().data()[1];

    auto fock_eig = array_ops::array_to_eigen(fock_recip);
    auto overlap_eig = array_ops::array_to_eigen(overlap);
    for (auto k = 0; k < k_size; ++k) {
        // Compute X = S^(-1/2)
        auto S = overlap_eig.block(0, k*tr0.extent(), tr0.extent(), tr0.extent());
        auto X = S.pow(-0.5);

        // Symmetrize Fock
        auto F = fock_eig.block(0, k*tr0.extent(), tr0.extent(), tr0.extent());
        F = (F + F.transpose().conjugate()) / 2.0;
        // F' = Xt.F.X
        Matrixc Xt = X.transpose().conjugate();
        auto XtF = Xt * F;
        auto Ft = XtF * X;

        // Diagonalize F'
        Eig::ComplexEigenSolver<Matrixc> comp_eig_solver(Ft);
        eps[k] = comp_eig_solver.eigenvalues();
        auto Ctemp = comp_eig_solver.eigenvectors();
        C[k] = X * Ctemp;
        // Sort eigenvalues and eigenvectors in ascending order
        sort_eigen(eps[k], C[k]);
    }

    Matrixc result_eig(tr0.extent(), tr1.extent());
    result_eig.setZero();
    auto R_size = 1 + idx_lattice(RD_max_(0), RD_max_(1), RD_max_(2), RD_max_);
    for (auto R = 0; R < R_size; ++R) {
        auto vec_R = R_vector(R, RD_max_);
        for (auto k = 0; k < k_size; ++k) {
            auto vec_k = k_vector(k);
            auto C_occ = C[k].leftCols(ndocc);
            auto D_real = C_occ.conjugate() * C_occ.transpose();
            auto exponent =
                    std::exp(I * vec_k.dot(vec_R)) / double(nk_(0) * nk_(1) * nk_(2));
            auto D_comp = exponent * D_real;
            result_eig.block(0, R*tr0.extent(), tr0.extent(), tr0.extent()) += D_comp;
        }
    }
    result = array_ops::eigen_to_array<Tile>(world_, result_eig, tr0, tr1);

    density_ = result;

    return result;
}


template <typename Tile, typename Policy>
void PeriodicAtomicIntegral<Tile, Policy>::sort_eigen(Vectorc &eigVal,
                                                      Matrixc &eigVec) {
  auto val = eigVal.real();

  // Sort by ascending eigenvalues
  std::vector<std::pair<double, int>> sortedVal;
  sortedVal.reserve(val.size());
  for (auto i = 0; i != val.size(); ++i) {
    auto pair = std::make_pair(val(i), i);
    sortedVal.push_back(pair);
  };
  std::sort(sortedVal.begin(), sortedVal.end());

  // Build sorted eigenvalues and eigenvectors
  Vectorc sortedEigVal(eigVal);
  Matrixc sortedEigVec(eigVec);
  for (auto i = 0; i != val.size(); ++i) {
    sortedEigVal(i) = eigVal(sortedVal[i].second);
    sortedEigVec.col(i) = eigVec.col(sortedVal[i].second);
  }

  eigVal = sortedEigVal;
  eigVec = sortedEigVec;
}

}  // mpqc namespace
}  // integrals namespace
#endif  // MPQC_PERIODIC_ATOMIC_INTEGRAL_H
