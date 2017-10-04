#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_

// set to 1, must have libint-2.4.0-beta.2
#define PRODUCE_PNO_MOLDEN_FILES 0
#if PRODUCE_PNO_MOLDEN_FILES
#include "libint2/lcao/molden.h"
#endif

#include "mpqc/chemistry/molecule/common.h"
#include "mpqc/chemistry/qc/cc/solvers.h"
#include "mpqc/chemistry/qc/lcao/factory/factory.h"
#include "mpqc/math/linalg/diagonal_array.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace lcao {
namespace cc {

namespace detail {
inline void print_local(TiledArray::World& world,
                        const std::shared_ptr<TiledArray::Pmap>& pmap) {
  for (ProcessID r = 0; r < world.size(); ++r) {
    world.gop.fence();
    if (r == world.rank()) {
      std::cout << r << ": { ";
      for (TiledArray::Pmap::const_iterator it = pmap->begin();
           it != pmap->end(); ++it)
        std::cout << *it << " ";
      std::cout << "}\n";
    }
  }
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> jacobi_update_t2_abij(
    const TA::DistArray<Tile, Policy>& r2_abij,
    const EigenVector<typename Tile::numeric_type>& ens_occ,
    const EigenVector<typename Tile::numeric_type>& ens_uocc) {
  auto denom = [&ens_occ, &ens_uocc](Tile& result_tile, const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];
    const auto i0 = result_tile.range().lobound()[2];
    const auto in = result_tile.range().upbound()[2];
    const auto j0 = result_tile.range().lobound()[3];
    const auto jn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::scalar_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens_uocc[a];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens_uocc[b];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens_occ[i];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens_occ[j];
            const auto e_iajb = e_i + e_j - e_a - e_b;
            const auto old = arg_tile[tile_idx];
            const auto result_abij = old / (e_iajb);
            const auto abs_result_abij = std::abs(result_abij);
            norm += abs_result_abij * abs_result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto delta_t2_abij = TA::foreach (r2_abij, denom);
  delta_t2_abij.world().gop.fence();
  return delta_t2_abij;
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> jacobi_update_t1_ai(
    const TA::DistArray<Tile, Policy>& r1_ai,
    const EigenVector<typename Tile::numeric_type>& ens_occ,
    const EigenVector<typename Tile::numeric_type>& ens_uocc) {
  auto denom = [&ens_occ, &ens_uocc](Tile& result_tile, const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto i0 = result_tile.range().lobound()[1];
    const auto in = result_tile.range().upbound()[1];

    auto tile_idx = 0;
    typename Tile::scalar_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens_uocc[a];
      for (auto i = i0; i < in; ++i, ++tile_idx) {
        const auto e_i = ens_occ[i];
        const auto e_ia = e_i - e_a;
        const auto old = arg_tile[tile_idx];
        const auto result_ai = old / (e_ia);
        const auto abs_result_ai = std::abs(result_ai);
        norm += abs_result_ai * abs_result_ai;
        result_tile[tile_idx] = result_ai;
      }
    }
    return std::sqrt(norm);
  };

  auto delta_t1_ai = TA::foreach (r1_ai, denom);
  delta_t1_ai.world().gop.fence();
  return delta_t1_ai;
}

/**
 *
 * @param r2_abij  T2 like Array with dimention "a,b,i,j", it should blocked by
 * V in a,b and 1 in i,j
 * @param F_occ_act  active occupied Fock matrix in canonical basis
 * @param F_pno_diag vector of diagonal Fock matrix in PNO basis, only local i,j
 * will be initialized
 * @param pnos vector of PNOs, only local i,j will be initialized
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> pno_jacobi_update_t2(
    const TA::DistArray<Tile, Policy>& r2_abij,
    const RowMatrix<typename Tile::numeric_type>& F_occ_act,
    const std::vector<EigenVector<typename Tile::numeric_type>>& F_pno_diag,
    const std::vector<RowMatrix<typename Tile::numeric_type>>& pnos) {
  auto update2 = [&F_occ_act, &F_pno_diag, &pnos](Tile& result_tile,
                                                  const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // determine i and j indices
    const auto i = arg_tile.range().lobound()[2];
    const auto j = arg_tile.range().lobound()[3];

    // Select appropriate matrix of PNOs
    auto ij = i * F_occ_act.cols() + j;
    Eigen::MatrixXd pno_ij = pnos[ij];

    // Extent data of tile
    const auto ext = arg_tile.range().extent_data();

    // Convert data in tile to Eigen::Map and transform to PNO basis
    const Eigen::MatrixXd r2_pno =
        pno_ij.transpose() *
        TA::eigen_map(arg_tile, ext[0] * ext[2], ext[1] * ext[3]) * pno_ij;

    // Create a matrix delta_t2_pno to hold updated values of delta_t2 in PNO
    // basis this matrix will then be back transformed to full basis before
    // being converted to a tile
    Eigen::MatrixXd delta_t2_pno = r2_pno;

    // Select correct vector containing diagonal elements of Fock matrix in
    // PNO basis
    const Eigen::VectorXd& ens_uocc = F_pno_diag[ij];

    // Determine number of PNOs
    const auto npno = ens_uocc.rows();

    // Determine number of uocc
    const auto nuocc = pno_ij.rows();

    // Select e_i and e_j
    const auto e_i = F_occ_act(i, i);
    const auto e_j = F_occ_act(j, j);

    for (auto a = 0; a < npno; ++a) {
      const auto e_a = ens_uocc[a];
      for (auto b = 0; b < npno; ++b) {
        const auto e_b = ens_uocc[b];
        const auto e_abij = e_i + e_j - e_a - e_b;
        const auto r_abij = r2_pno(a, b);
        delta_t2_pno(a, b) = r_abij / e_abij;
      }
    }

    // Back transform delta_t2_pno to full space
    Eigen::MatrixXd delta_t2_full = pno_ij * delta_t2_pno * pno_ij.transpose();

    // Convert delta_t2_full to tile and compute norm
    typename Tile::scalar_type norm = 0.0;
    for (auto r = 0; r < nuocc; ++r) {
      for (auto c = 0; c < nuocc; ++c) {
        const auto idx = r * nuocc + c;
        const auto elem = delta_t2_full(r, c);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[idx] = elem;
      }
    }

    return std::sqrt(norm);
  };
  auto delta_t2_abij = TA::foreach (r2_abij, update2);
  delta_t2_abij.world().gop.fence();
  return delta_t2_abij;
}

/**
 *
 * @param r1_ai  T1 like Array with dimention "a,i", it should blocked by V in
 * a,b and 1 in i,j
 * @param F_occ_act  active occupied Fock matrix in canonical basis
 * @param F_pno_diag vector of diagonal Fock matrix in PNO basis, only local i,j
 * will be initialized
 * @param osvs vector of OSVs, only local i will be initialized
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> pno_jacobi_update_t1(
    const TA::DistArray<Tile, Policy>& r1_ai,
    const RowMatrix<typename Tile::numeric_type>& F_occ_act,
    const std::vector<EigenVector<typename Tile::numeric_type>>& F_osv_diag,
    const std::vector<RowMatrix<typename Tile::numeric_type>>& osvs) {
  auto update1 = [&F_occ_act, &F_osv_diag, &osvs](Tile& result_tile,
                                                  const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // determine i index
    const auto i = arg_tile.range().lobound()[1];

    // Select appropriate matrix of OSVs
    Eigen::MatrixXd osv_i = osvs[i];

    // Extent data of tile
    const auto ext = arg_tile.range().extent_data();

    // Convert data in tile to Eigen::Map and transform to OSV basis
    const Eigen::VectorXd r1_osv =
        osv_i.transpose() * TA::eigen_map(arg_tile, ext[0], ext[1]);

    // Create a matrix delta_t1_osv to hold updated values of delta t1 in OSV
    // basis. This matrix will then be back transformed to full basis before
    // being converted to a tile
    Eigen::VectorXd delta_t1_osv = r1_osv;

    // Select correct vector containing diagonal elements of Fock matrix in
    // OSV basis
    const Eigen::VectorXd& ens_uocc = F_osv_diag[i];

    // Determine number of OSVs
    //      const auto nosv = ens_uocc.rows();
    const auto nosv = osv_i.cols();

    // Determine number of uocc
    const auto nuocc = osv_i.rows();

    // Select e_i
    const auto e_i = F_occ_act(i, i);

    for (auto a = 0; a < nosv; ++a) {
      const auto e_a = ens_uocc[a];
      const auto e_ai = e_i - e_a;
      const auto r_ai = r1_osv(a);
      delta_t1_osv(a) = r_ai / e_ai;
    }

    // Back transform delta_t1_osv to full space
    Eigen::VectorXd delta_t1_full = osv_i * delta_t1_osv;

    // Convert delta_t1_full to tile and compute norm
    typename Tile::scalar_type norm = 0.0;
    for (auto r = 0; r < nuocc; ++r) {
      const auto elem = delta_t1_full(r);
      const auto abs_elem = std::abs(elem);
      norm += abs_elem * abs_elem;
      result_tile[r] = elem;
    }

    return std::sqrt(norm);
  };

  auto delta_t1_ai = TA::foreach (r1_ai, update1);
  delta_t1_ai.world().gop.fence();
  return delta_t1_ai;
}

/**
 *
 * @param abij T2 like Array with dimention "a,b,i,j", it should blocked by V in
 * a,b and 1 in i,j
 * @param pnos vector of PNOs, only local i,j will be initialized
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> pno_transform_abij(
    const TA::DistArray<Tile, Policy>& abij,
    const std::vector<RowMatrix<typename Tile::numeric_type>>& pnos) {
  std::size_t nocc_act = abij.trange().dim(2).extent();

  auto tform = [&pnos, nocc_act](Tile& result_tile, const Tile& arg_tile) {

    // determine i and j indices
    const auto i = arg_tile.range().lobound()[2];
    const auto j = arg_tile.range().lobound()[3];

    // Select appropriate matrix of PNOs
    const auto ij = i * nocc_act + j;
    Eigen::MatrixXd pno_ij = pnos[ij];
    const auto nuocc = pno_ij.rows();
    const auto npno = pno_ij.cols();

    // Convert data in tile to Eigen::Map and transform to PNO basis
    const Eigen::MatrixXd result_eig =
        pno_ij.transpose() * TA::eigen_map(arg_tile, nuocc, nuocc) * pno_ij;

    // Convert result_eig to tile and compute norm
    result_tile = Tile(TA::Range{npno, npno, 1l, 1l});
    typename Tile::scalar_type norm = 0.0;
    for (auto r = 0; r < npno; ++r) {
      for (auto c = 0; c < npno; ++c) {
        const auto idx = r * npno + c;
        const auto elem = result_eig(r, c);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[idx] = elem;
      }
    }

    return std::sqrt(norm);
  };

  auto result = TA::foreach (abij, tform);
  result.world().gop.fence();
  return result;
}

/**
 *
 * @param ai  T1 like Array with dimention "a,i", it should blocked by V in a,b
 * and 1 in i,j
 * @param osvs vector of OSVs, only local i will be initialized
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> osv_transform_ai(
    const TA::DistArray<Tile, Policy>& ai,
    const std::vector<RowMatrix<typename Tile::numeric_type>>& osvs) {
  auto tform = [&osvs](Tile& result_tile, const Tile& arg_tile) {

    // determine i index
    const auto i = arg_tile.range().lobound()[1];

    // Select appropriate matrix of OSVs
    Eigen::MatrixXd osv_i = osvs[i];
    const auto nuocc = osv_i.rows();
    const auto nosv = osv_i.cols();

    // Convert data in tile to Eigen::Map and transform to OSV basis
    const Eigen::MatrixXd result_eig =
        osv_i.transpose() * TA::eigen_map(arg_tile, nuocc, 1);

    // Convert result_eig to tile and compute norm
    result_tile = Tile(TA::Range{nosv, 1l});
    typename Tile::scalar_type norm = 0.0;
    for (auto r = 0; r < nosv; ++r) {
      const auto elem = result_eig(r, 0);
      const auto abs_elem = std::abs(elem);
      norm += abs_elem * abs_elem;
      result_tile[r] = elem;
    }

    return std::sqrt(norm);
  };

  auto result = TA::foreach (ai, tform);
  result.world().gop.fence();
  return result;
}

/**
 *
 * @param t2   T2 like array with dimention "a,b,i,j"
 * @param i_block array "i,j" reblock occupied dimention i to j
 * @param a_block array "a,b" that reblock unoccupied dimention a to b
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> reblock_t2(
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& i_block,
    const TA::DistArray<Tile, Policy>& a_block) {
  // Reblock T2
  TA::DistArray<Tile, Policy> result;
  result("an,bn,in,jn") = t2("a,b,i,j") * i_block("j,jn") * i_block("i,in") *
                          a_block("b,bn") * a_block("a,an");

  return result;
}

/**
 *
 * @param t2   T2 like array with dimention "a,b,i,j"
 * @param i_block array "i,j" reblock occupied dimention j to i
 * @param a_block array "a,b" that reblock unoccupied dimention b to a
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> unblock_t2(
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& i_block,
    const TA::DistArray<Tile, Policy>& a_block) {
  // Reblock T2
  TA::DistArray<Tile, Policy> result;
  result("an,bn,in,jn") = t2("a,b,i,j") * i_block("jn,j") * i_block("in,i") *
                          a_block("bn,b") * a_block("an,a");

  return result;
}

/**
 *
 * @param t1  T1 like array with dimention "a,i"
 * @param i_block array "i,j" reblock occupied dimention i to j
 * @param a_block array "a,b" that reblock unoccupied dimention a to b
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> reblock_t1(
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& i_block,
    const TA::DistArray<Tile, Policy>& a_block) {
  // Create target_pmap to match OSVs pma
  auto nocc_act = i_block.range().extent_data()[1];
  auto target_pmap = Policy::default_pmap(t1.world(), nocc_act);

  // Reblock T1
  TA::DistArray<Tile, Policy> result;
  result("an,in") = (t1("a,i") * i_block("i,in") * a_block("a,an")).set_pmap(target_pmap);
  return result;
}

/**
 *
 * @param t1  T1 like array with dimention "a,i"
 * @param i_block array "i,j" reblock occupied dimention j to i
 * @param a_block array "a,b" that reblock unoccupied dimention b to a
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> unblock_t1(
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& i_block,
    const TA::DistArray<Tile, Policy>& a_block) {
  // Reblock T1
  TA::DistArray<Tile, Policy> result;
  result("an,in") = t1("a,i") * i_block("in,i") * a_block("an,a");
  return result;
}

/**
 *
 * @param t2
 * @param pnos
 * @param F_pno_diags
 * @param osvs
 * @param F_osv_diags
 */
template <typename Tile, typename Policy>
void construct_pno(
    const TA::DistArray<Tile, Policy>& t2,
    const RowMatrix<typename Tile::numeric_type>& F_uocc,
    double tpno,
    double tosv,
    std::vector<RowMatrix<typename Tile::numeric_type>>& pnos,
    std::vector<int>& npnos,
    std::vector<EigenVector<typename Tile::numeric_type>>& F_pno_diag,
    std::vector<RowMatrix<typename Tile::numeric_type>>& osvs,
    std::vector<int>& nosvs,
    std::vector<EigenVector<typename Tile::numeric_type>>& F_osv_diag,
    bool pno_canonical = false) {
  using Matrix = RowMatrix<typename Tile::numeric_type>;
  using Vector = EigenVector<typename Tile::numeric_type>;

  auto& world = t2.world();
  std::size_t nocc_act = t2.trange().dim(2).extent();
  std::size_t nuocc = t2.trange().dim(0).extent();

  // For storing PNOs and and the Fock matrix in the PNO basis
  npnos.resize(nocc_act * nocc_act, 0);
  pnos.resize(nocc_act * nocc_act);
  F_pno_diag.resize(nocc_act * nocc_act);


  // For storing OSVs (PNOs when i = j) and the Fock matrix in
  // the OSV basis
  nosvs.resize(nocc_act, 0);
  osvs.resize(nocc_act);
  F_osv_diag.resize(nocc_act);

  /// Step (2): Transform T to D

  // lambda function to transform T to D; implement using a for_each
  auto form_D = [nuocc](Tile& result_tile, const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // Get values of i,j and compute delta_ij
    const int i = arg_tile.range().lobound()[2];
    const int j = arg_tile.range().lobound()[3];

    auto delta_ij = (i == j) ? 1 : 0;

    // Form T_ij matrix from arg_tile
    Matrix T_ij = TA::eigen_map(arg_tile, nuocc, nuocc);

    // Form T_ji matrix from T_ij
    Matrix T_ji = T_ij.transpose();

    // Form D_ij from T_ij and T_ji
    Matrix D_ij = (1.0 / (1 + delta_ij)) * (4 * T_ji * T_ij - 2 * T_ij * T_ij +
                                            4 * T_ij * T_ji - 2 * T_ji * T_ji);

    // Transform D_ij into a tile
    auto norm = 0.0;
    for (int a = 0, tile_idx = 0; a != nuocc; ++a) {
      for (int b = 0; b != nuocc; ++b, ++tile_idx) {
        const auto elem = D_ij(a, b);
        const auto abs_result = std::abs(elem);
        norm += (abs_result * abs_result);
        result_tile[tile_idx] = elem;
      }
    }
    return std::sqrt(norm);
  };  // form_D

  auto D = TA::foreach (t2, form_D);
  D.world().gop.fence();

  /// Step (3): Form PNO_ij from D_ij
  // Diagonalize each D_ij matrix to get the PNOs and occupation
  // numbers

  // Lambda function to form osvs
  auto form_OSV = [&D, &F_osv_diag,
                   &F_uocc, &nosvs, tosv, nuocc, nocc_act,
                   pno_canonical](TA::World& world){

    Eigen::SelfAdjointEigenSolver<Matrix> es;

    std::vector<Matrix> osvs_list(nocc_act);

    auto reblock_r1_pmap = Policy::default_pmap(world, nocc_act);

    for (ProcessID r = 0; r < world.size(); ++r) {
      std::vector<int> local_tiles;
      world.gop.fence();
      if (r == world.rank()) {
        for (TiledArray::Pmap::const_iterator it = reblock_r1_pmap->begin();
             it != reblock_r1_pmap->end(); ++it){
          local_tiles.push_back(*it);
        }

        for (int idx = 0; idx != local_tiles.size(); ++idx) {
          auto i = local_tiles[idx];

          Tile D_ii = D.find({0,0,i,i}).get();
          Matrix D_ii_mat = TA::eigen_map(D_ii, nuocc, nuocc);


          // Diagonalize D_ii_mat
          es.compute(D_ii_mat);
          Matrix osv_ii = es.eigenvectors();
          auto occ_ii = es.eigenvalues();

          // Truncate osvs_list
          size_t osvdrop = 0;
          if (tosv != 0.0) {
            for (size_t k = 0; k != occ_ii.rows(); ++k) {
              if (!(occ_ii(k) >= tosv))
                ++osvdrop;
              else
                break;
            }
          }
          const auto nosv = nuocc - osvdrop;
          nosvs[i] = nosv;

          if (nosv == 0) {  // all OSV truncated indicates total nonsense
            throw LimitExceeded<size_t>("all osvs_list truncated", __FILE__, __LINE__, 1,
                                        0);
          }

          // Store truncated OSV mat in osvs_list
          Matrix osv_trunc = osv_ii.block(0, osvdrop, nuocc, nosv);
          osvs_list[i] = osv_trunc;

          // Transform F to OSV space
          Matrix F_osv_i = osv_trunc.transpose() * F_uocc * osv_trunc;

          // Store just the diagonal elements of F_osv_i
          F_osv_diag[i] = F_osv_i.diagonal();

          /////// Transform osvs_list to canonical osvs_list if pno_canonical_ == true
          if (pno_canonical) {
            // Compute eigenvectors of F in OSV space
            es.compute(F_osv_i);
            Matrix osv_transform_i = es.eigenvectors();

            // Transform osv_i to canonical OSV space: osv_i -> can_osv_i
            Matrix can_osv_i = osv_trunc * osv_transform_i;

            // Replace standard with canonical osvs_list
            osvs_list[i] = can_osv_i;
            F_osv_diag[i] = es.eigenvalues();
          }  // pno_canonical
        }
      }
    }
    return osvs_list;
  };

  osvs = form_OSV(world);
  // The following fence is unnecessary IF OSVs are only used locally
  // world.gop.fence();

  // Lambda function to form PNOs; implement using a for_each
  auto form_PNO = [&pnos, &F_pno_diag, &osvs, &F_osv_diag, &F_uocc, &npnos,
                   &nosvs, tpno, tosv, nuocc, nocc_act,
                   pno_canonical](Tile& result_tile, const Tile& arg_tile) {

    Eigen::SelfAdjointEigenSolver<Matrix> es;

    result_tile = Tile(arg_tile.range());

    // Get values of i,j and compute ij
    const int i = arg_tile.range().lobound()[2];
    const int j = arg_tile.range().lobound()[3];
    const int ij = i * nocc_act + j;

    // Form D_ij matrix from arg_tile
    Matrix D_ij = TA::eigen_map(arg_tile, nuocc, nuocc);

    // Diagonalize D_ij
    es.compute(D_ij);
    Matrix pno_ij = es.eigenvectors();
    auto occ_ij = es.eigenvalues();

    // Determine number of PNOs to be dropped
    std::size_t pnodrop = 0;

    if (tpno != 0.0) {
      for (std::size_t k = 0; k != occ_ij.rows(); ++k) {
        if (!(occ_ij(k) >= tpno))
          ++pnodrop;
        else
          break;
      }
    }

    // Calculate the number of PNOs kept and store in vector npnos
    // for calculating average npno later
    const auto npno = nuocc - pnodrop;
    npnos[ij] = npno;

    // Truncate PNO matrix and store in pnos_

    // Declare matrix pno_trunc to hold truncated set of PNOs
    Matrix pno_trunc;

    // If npno = 0, substitute a single fake "PNO" with all coefficients
    // equal to zero. All other code will behave the same way
    if (npno == 0) {
      // resize pno_trunc to be size nocc_act x 1
      pno_trunc.resize(nuocc, 1);

      // Create a zero matrix of size nuocc x 1
      Matrix pno_zero = Matrix::Zero(nuocc, 1);

      // Set pno_trunc eqaul to pno_zero matrix
      pno_trunc = pno_zero;
    }

    // If npno != zero, use actual zet of truncated PNOs
    else {
      pno_trunc.resize(nuocc, npno);
      pno_trunc = pno_ij.block(0, pnodrop, nuocc, npno);
    }

    // Store truncated PNOs
    pnos[ij] = pno_trunc;

    // Transform F to PNO space
    Matrix F_pno_ij = pno_trunc.transpose() * F_uocc * pno_trunc;

    // Store just the diagonal elements of F_pno_ij
    F_pno_diag[ij] = F_pno_ij.diagonal();

    ///// Transform PNOs to canonical PNOs if pno_canonical_ == true

    if (pno_canonical && npno > 0) {
      //          if (pno_canonical_ == "true") {
      // Compute eigenvectors of F in PNO space
      es.compute(F_pno_ij);
      Matrix pno_transform_ij = es.eigenvectors();

      // Transform pno_ij to canonical PNO space; pno_ij -> can_pno_ij
      Matrix can_pno_ij = pno_trunc * pno_transform_ij;

      // Replace standard with canonical PNOs
      pnos[ij] = can_pno_ij;
      F_pno_diag[ij] = es.eigenvalues();
    }  // pno_canonical

    // Transform D_ij into a tile
    auto norm = 0.0;
    for (int a = 0, tile_idx = 0; a != nuocc; ++a) {
      for (int b = 0; b != nuocc; ++b, ++tile_idx) {
        const auto elem = D_ij(a, b);
        const auto abs_result = std::abs(elem);
        norm += abs_result;
        result_tile[tile_idx] = elem;
      }
    }
    return std::sqrt(norm);
  };  // form_PNO

  auto D_prime = TA::foreach (D, form_PNO);
  world.gop.fence();

  // Sum together vectors of npnos and nosvs on each node
  world.gop.sum(npnos.data(), npnos.size());
  world.gop.sum(nosvs.data(), nosvs.size());

  // Compute and print average number of OSVs per pair
  if (D_prime.world().rank() == 0) {
    auto tot_osv = 0;
    for (int i = 0; i != nosvs.size(); ++i) {
      tot_osv += nosvs[i];
    }

    auto ave_nosv = tot_osv / nocc_act;
    ExEnv::out0() << "The average number of OSVs per pair is " << ave_nosv
                  << std::endl;

    // Compute and print average number of PNOs per pair

    auto tot_pno = 0;
    for (int i = 0; i != npnos.size(); ++i) {
      tot_pno += npnos[i];
    }

    auto ave_npno = tot_pno / (nocc_act * nocc_act);
    ExEnv::out0() << "The average number of PNOs per pair is " << ave_npno
                  << std::endl;
  }  // end if D_prime_.world == 0

};  // construct_pno


// Reconstruct PNOs using updated T2
template <typename Tile, typename Policy>
void reconstruct_pno(
    const TA::DistArray<Tile, Policy>& t2,
    const RowMatrix<typename Tile::numeric_type>& F_uocc,
    std::vector<int>& npnos,
    std::vector<int>& nosvs,
    std::vector<RowMatrix<typename Tile::numeric_type>>& pnos,
    std::vector<EigenVector<typename Tile::numeric_type>>& F_pno_diag,
    std::vector<RowMatrix<typename Tile::numeric_type>>& osvs,
    std::vector<EigenVector<typename Tile::numeric_type>>& F_osv_diag,
    bool pno_canonical = false) {
  using Matrix = RowMatrix<typename Tile::numeric_type>;
  using Vector = EigenVector<typename Tile::numeric_type>;

  auto& world = t2.world();
  std::size_t nocc_act = t2.trange().dim(2).extent();
  std::size_t nuocc = t2.trange().dim(0).extent();

  // For storing PNOs and and the Fock matrix in the PNO basis
  pnos.resize(nocc_act * nocc_act);
  F_pno_diag.resize(nocc_act * nocc_act);


  // For storing OSVs (PNOs when i = j) and the Fock matrix in
  // the OSV basis
  osvs.resize(nocc_act);
  F_osv_diag.resize(nocc_act);

  /// Step (2): Transform T to D

  // lambda function to transform T to D; implement using a for_each
  auto form_D = [nuocc](Tile& result_tile, const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // Get values of i,j and compute delta_ij
    const int i = arg_tile.range().lobound()[2];
    const int j = arg_tile.range().lobound()[3];

    auto delta_ij = (i == j) ? 1 : 0;

    // Form T_ij matrix from arg_tile
    Matrix T_ij = TA::eigen_map(arg_tile, nuocc, nuocc);

    // Form T_ji matrix from T_ij
    Matrix T_ji = T_ij.transpose();

    // Form D_ij from T_ij and T_ji
    Matrix D_ij = (1.0 / (1 + delta_ij)) * (4 * T_ji * T_ij - 2 * T_ij * T_ij +
                                            4 * T_ij * T_ji - 2 * T_ji * T_ji);

    // Transform D_ij into a tile
    auto norm = 0.0;
    for (int a = 0, tile_idx = 0; a != nuocc; ++a) {
      for (int b = 0; b != nuocc; ++b, ++tile_idx) {
        const auto elem = D_ij(a, b);
        const auto abs_result = std::abs(elem);
        norm += (abs_result * abs_result);
        result_tile[tile_idx] = elem;
      }
    }
    return std::sqrt(norm);
  };  // form_D

  auto D = TA::foreach (t2, form_D);
  D.world().gop.fence();


  /// Step (3): Form PNO_ij from D_ij
  // Diagonalize each D_ij matrix to get the PNOs and occupation
  // numbers

  // Lambda function to form osvs
  auto form_OSV = [&D, &F_osv_diag,
                   &F_uocc, &nosvs, nuocc, nocc_act,
                   pno_canonical](TA::World& world){

    Eigen::SelfAdjointEigenSolver<Matrix> es;

    std::vector<Matrix> osvs_list(nocc_act);

    auto reblock_r1_pmap = Policy::default_pmap(world, nocc_act);

    for (ProcessID r = 0; r < world.size(); ++r) {
      std::vector<int> local_tiles;
      world.gop.fence();
      if (r == world.rank()) {
        for (TiledArray::Pmap::const_iterator it = reblock_r1_pmap->begin();
             it != reblock_r1_pmap->end(); ++it){
          local_tiles.push_back(*it);
        }

        for (int idx = 0; idx != local_tiles.size(); ++idx) {
          auto i = local_tiles[idx];

          Tile D_ii = D.find({0,0,i,i}).get();
          Matrix D_ii_mat = TA::eigen_map(D_ii, nuocc, nuocc);


          // Diagonalize D_ii_mat
          es.compute(D_ii_mat);
          Matrix osv_ii = es.eigenvectors();

          // Use the same number of OSVs per i as used when originally forming OSVs
          auto nosv = nosvs[i];
          auto osvdrop = nuocc - nosv;

          if (nosv == 0) {  // all OSV truncated indicates total nonsense
            throw LimitExceeded<size_t>("all osvs_list truncated", __FILE__, __LINE__, 1,
                                        0);
          }

          // Store truncated OSV mat in osvs_list
          Matrix osv_trunc = osv_ii.block(0, osvdrop, nuocc, nosv);
          osvs_list[i] = osv_trunc;

          // Transform F to OSV space
          Matrix F_osv_i = osv_trunc.transpose() * F_uocc * osv_trunc;

          // Store just the diagonal elements of F_osv_i
          F_osv_diag[i] = F_osv_i.diagonal();

          /////// Transform osvs_list to canonical osvs_list if pno_canonical_ == true
          if (pno_canonical) {
            // Compute eigenvectors of F in OSV space
            es.compute(F_osv_i);
            Matrix osv_transform_i = es.eigenvectors();

            // Transform osv_i to canonical OSV space: osv_i -> can_osv_i
            Matrix can_osv_i = osv_trunc * osv_transform_i;

            // Replace standard with canonical osvs_list
            osvs_list[i] = can_osv_i;
            F_osv_diag[i] = es.eigenvalues();
          }  // pno_canonical
        }
      }
    }
    return osvs_list;
  };

  osvs = form_OSV(world);
  // The following fence is unnecessary IF OSVs are only used locally
  // world.gop.fence();

  // Lambda function to form PNOs; implement using a for_each
  auto form_PNO = [&pnos, &F_pno_diag, &osvs, &F_osv_diag, &F_uocc, &npnos,
                   &nosvs, nuocc, nocc_act,
                   pno_canonical](Tile& result_tile, const Tile& arg_tile) {

    Eigen::SelfAdjointEigenSolver<Matrix> es;

    result_tile = Tile(arg_tile.range());

    // Get values of i,j and compute ij
    const int i = arg_tile.range().lobound()[2];
    const int j = arg_tile.range().lobound()[3];
    const int ij = i * nocc_act + j;

    // Form D_ij matrix from arg_tile
    Matrix D_ij = TA::eigen_map(arg_tile, nuocc, nuocc);


    // Diagonalize D_ij
    es.compute(D_ij);
    Matrix pno_ij = es.eigenvectors();
    auto occ_ij = es.eigenvalues();


    // Use the same number of PNOs were i,j pair as when first forming PNOs
    auto npno = npnos[ij];
    auto pnodrop = nuocc - npno;

    // Truncate PNO matrix and store in pnos_

    // Declare matrix pno_trunc to hold truncated set of PNOs
    Matrix pno_trunc;

    // If npno = 0, substitute a single fake "PNO" with all coefficients
    // equal to zero. All other code will behave the same way
    if (npno == 0) {
      // resize pno_trunc to be size nocc_act x 1
      pno_trunc.resize(nuocc, 1);

      // Create a zero matrix of size nuocc x 1
      Matrix pno_zero = Matrix::Zero(nuocc, 1);

      // Set pno_trunc eqaul to pno_zero matrix
      pno_trunc = pno_zero;
    }

    // If npno != zero, use actual zet of truncated PNOs
    else {
      pno_trunc.resize(nuocc, npno);
      pno_trunc = pno_ij.block(0, pnodrop, nuocc, npno);
    }

    // Store truncated PNOs
    pnos[ij] = pno_trunc;

    // Transform F to PNO space
    Matrix F_pno_ij = pno_trunc.transpose() * F_uocc * pno_trunc;

    // Store just the diagonal elements of F_pno_ij
    F_pno_diag[ij] = F_pno_ij.diagonal();

    ///// Transform PNOs to canonical PNOs if pno_canonical_ == true

    if (pno_canonical && npno > 0) {
      //          if (pno_canonical_ == "true") {
      // Compute eigenvectors of F in PNO space
      es.compute(F_pno_ij);
      Matrix pno_transform_ij = es.eigenvectors();

      // Transform pno_ij to canonical PNO space; pno_ij -> can_pno_ij
      Matrix can_pno_ij = pno_trunc * pno_transform_ij;

      // Replace standard with canonical PNOs
      pnos[ij] = can_pno_ij;
      F_pno_diag[ij] = es.eigenvalues();
    }  // pno_canonical

    // Transform D_ij into a tile
    auto norm = 0.0;
    for (int a = 0, tile_idx = 0; a != nuocc; ++a) {
      for (int b = 0; b != nuocc; ++b, ++tile_idx) {
        const auto elem = D_ij(a, b);
        const auto abs_result = std::abs(elem);
        norm += abs_result;
        result_tile[tile_idx] = elem;
      }
    }
    return std::sqrt(norm);
  };  // form_PNO

  auto D_prime = TA::foreach (D, form_PNO);
  world.gop.fence();

};  // reconstruct_pno



}  // namespace detail

// // JacobiDIISSolver updates the CC T amplitudes using standard Jacobi+DIIS
template <typename T>
class JacobiDIISSolver : public ::mpqc::cc::DIISSolver<T, T> {
 public:
  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object; it will be queried for all keywords of ::mpqc::cc::DIISSolver<T,T> .
   */
  // clang-format on

  JacobiDIISSolver(const KeyVal& kv,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> f_ii,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> f_aa)
      : ::mpqc::cc::DIISSolver<T, T>(kv) {
    std::swap(f_ii_, f_ii);
    std::swap(f_aa_, f_aa);
  }
  virtual ~JacobiDIISSolver() = default;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 1> f_ii_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> f_aa_;

  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    t1("a,i") += detail::jacobi_update_t1_ai(r1, f_ii_, f_aa_)("a,i");
    t2("a,b,i,j") += detail::jacobi_update_t2_abij(r2, f_ii_, f_aa_)("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }
};

/// PNOSolver updates the CC T amplitudes using standard Jacobi+DIIS in PNO
/// space
/// @warning This class assumes that the 1- and 2-body amplitudes/residuals
///          given to Solver::update() are laid out as "a,i" and "a,b,i,j",
///          respectively
template <typename T, typename DT>
class PNOSolver : public ::mpqc::cc::DIISSolver<T, T>,
                  public madness::WorldObject<PNOSolver<T, DT>> {
 public:
  using Tile = typename T::value_type;
  using Matrix = RowMatrix<typename Tile::numeric_type>;
  using Vector = EigenVector<typename Tile::numeric_type>;

  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object; it will be queried for all keywords of ::mpqc::cc::DIISSolver , as well
   * as the following additional keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | pno_method | string | standard | The PNO construction method. Valid values are: \c standard . |
   * | tpno | double | 1e-8 | The PNO construction threshold. This non-negative integer specifies the screening threshold for the eigenvalues of the pair density. Setting this to zero will cause the full (untruncated) set of PNOs to be used. |
   * | tosv | double | 1e-9 | The OSV construction threshold. This non-negative integer specifies the screening threshold for the eigenvalues of the pair density of the diagonal pairs. Setting this to zero will cause the full (untruncated) set of OSVs to be used. |
   * | pno_canonical | bool | false | Whether or not to canonicalize the PNOs and OSVs |
   * | tiling_method | string | flexible | How the basis set is tiled. Valid values are: \c flexible , \c rigid . |
   * | interval | int | 10 | Every nth iteration, PNOs are recomputed |
   * | residual_thresh | double | 1e-10 | Once the residual value is less than the threshold value, update_pno set to false |
   */
  // clang-format on
  PNOSolver(const KeyVal& kv, Factory<T, DT>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv),
        madness::WorldObject<PNOSolver<T, DT>>(factory.world()),
        factory_(factory),
        pno_method_(kv.value<std::string>("pno_method", "standard")),
        pno_canonical_(kv.value<bool>("pno_canonical", false)),
        update_pno_(kv.value<bool>("update_pno", false)),
        tiling_method_(kv.value<std::string>(
            "tiling_method", "flexible")),  // allows user to select whether
                                            // rigid or flexible tiling code is
                                            // employed
        tpno_(kv.value<double>("tpno", 1.e-8)),
        tosv_(kv.value<double>("tosv", 1.e-9)),
        interval_(kv.value<int>("interval", 10)),
        residual_thresh_(kv.value<double>("residual_thresh", 1e-10)){
    // part of WorldObject initialization
    this->process_pending();

    // compute and store PNOs truncated with threshold tpno_
    // store PNOs for diagonal pair as OSVs truncated with threshold tosv_

    // Determine value of tiling_method_
    std::string tiling_method_ = (kv.value<std::string>("tiling_method"));

    // If tiling_method_ = rigid, employ original PNO code
    if (tiling_method_ == "rigid") {
      // Check that tiling is done appropriately
      if (kv.exists("occ_block_size")) {
        int occ_block_size_ = (kv.value<int>("occ_block_size"));
        if (occ_block_size_ != 1) {
          throw InputError(
              "occ_block_size must be set to 1 in the input file.");
        }
      } else {
        throw InputError("occ_block_size was not specified in the input file.");
      }

      if (kv.exists("unocc_block_size")) {
        int unocc_block_size_ = (kv.value<int>("unocc_block_size"));
        if (unocc_block_size_ < 1000000000) {
          throw InputError(
              "unocc_block_size must be greater than or equal to 1000000000 in "
              "the input file.");
        }
      } else {
        throw InputError(
            "unocc_block_size was not specified in the input file.");
      }

      /// !!! Original PNO formation code !!! ///
      // Do not delete! //

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

      auto& fac = factory_;
      auto& world = fac.world();
      auto& ofac = fac.orbital_registry();

      auto nocc = ofac.retrieve("m").rank();
      auto nocc_act = ofac.retrieve("i").rank();
      nocc_act_ = nocc_act;
      auto nuocc = ofac.retrieve("a").rank();
      nuocc_ = nuocc;
      auto nfzc = nocc - nocc_act;

      // Form Fock array
      auto F = fac.compute(L"<p|F|q>[df]");

      // Select just diagonal elements of Fock aray and transform
      // to Eigen vector; use for computing PNOs
      Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
      auto eps_o = eps_p.segment(nfzc, nocc_act);
      auto eps_v = eps_p.tail(nuocc);

      // Transform entire Fock array to Eigen Matrix
      Eigen::MatrixXd F_all = array_ops::array_to_eigen(F);

      // Select just the occupied portion of the Fock matrix
      F_occ_act_ = F_all.block(nfzc, nfzc, nocc_act, nocc_act);

      // Select just the unoccupied portion of the Fock matrix
      Eigen::MatrixXd F_uocc = F_all.block(nocc, nocc, nuocc, nuocc);

      // Compute all K_aibj
      // auto K = fac.compute(L"(a b|G|i j)");
      auto K = fac.compute(L"<a b|G|i j>[df]");
      const auto ktrange = K.trange();

      // zero out amplitudes
      if (!T_.is_initialized()) {
        T_ = T(world, K.trange(), K.shape());
        T_.fill(0.0);
      }

      // For storing PNOs and and the Fock matrix in the PNO basis
      pnos_.resize(nocc_act * nocc_act);
      F_pno_diag_.resize(nocc_act * nocc_act);

      // For storing OSVs (PNOs when i = j) and the Fock matrix in
      // the OSV basis
      osvs_.resize(nocc_act);
      F_osv_diag_.resize(nocc_act);

#if PRODUCE_PNO_MOLDEN_FILES
      // prepare to Molden
      const auto libint2_atoms = to_libint_atom(fac.atoms()->atoms());
      const auto C_i_eig =
          array_ops::array_to_eigen(ofac.retrieve("i").coefs());
      const auto C_a_eig =
          array_ops::array_to_eigen(ofac.retrieve("a").coefs());
      const auto libint2_shells =
          fac.basis_registry()->retrieve(L"μ")->flattened_shells();

      // write out active occupied orbitals
      auto occs = Eigen::VectorXd::Constant(C_i_eig.cols(), 2.0);
      libint2::molden::Export xport(libint2_atoms, libint2_shells, C_i_eig,
                                    occs);
      xport.write("occ.molden");
#endif

      // Keep running sum of all PNOs in systems for computing average
      auto total_pno_sum = 0;

      // Loop over each pair of occupieds to form PNOs
      for (int i = 0; i < nocc_act; ++i) {
        double eps_i = eps_o[i];

        for (int j = 0; j < nocc_act; ++j) {
          double eps_j = eps_o[j];
          int delta_ij = (i == j) ? 1 : 0;
          std::array<int, 4> tile_ij = {{0, 0, i, j}};
          std::array<int, 4> tile_ji = {{0, 0, j, i}};
          const auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
          const auto ord_ji = ktrange.tiles_range().ordinal(tile_ji);
          TA::TensorD K_ij = K.find(ord_ij);
          TA::TensorD K_ji = K.find(ord_ji);
          auto ext_ij = K_ij.range().extent_data();
          auto ext_ji = K_ji.range().extent_data();
          Eigen::MatrixXd K_ij_mat =
              TA::eigen_map(K_ij, ext_ij[0] * ext_ij[2], ext_ij[1] * ext_ij[3]);
          Eigen::MatrixXd K_ji_mat =
              TA::eigen_map(K_ji, ext_ji[0] * ext_ji[2], ext_ji[1] * ext_ji[3]);

          Eigen::MatrixXd T_ij(nuocc, nuocc);
          Eigen::MatrixXd T_ji(nuocc, nuocc);
          Eigen::MatrixXd T_tilde_ij(nuocc, nuocc);

          for (int a = 0; a < nuocc; ++a) {
            double eps_a = eps_v[a];
            for (int b = 0; b < nuocc; ++b) {
              double eps_b = eps_v[b];

              T_ij(a, b) = -K_ij_mat(a, b) / (eps_a + eps_b - eps_i - eps_j);
              T_ji(a, b) = -K_ji_mat(a, b) / (eps_a + eps_b - eps_i - eps_j);
            }
          }

          // Eq. 23, JCP 128 034106 (2013)
          T_tilde_ij = 4 * T_ij - 2 * T_ji;
          Eigen::MatrixXd D_ij =
              (T_tilde_ij.transpose() * T_ij + T_tilde_ij * T_ij.transpose()) /
              (1.0 + delta_ij);

          // Diagonalize D_ij to get PNOs and corresponding occupation numbers.
          es.compute(D_ij);
          Eigen::MatrixXd pno_ij = es.eigenvectors();
          auto occ_ij = es.eigenvalues();

          //         std::cout << "D_" << i << j << " mat:\n" << D_ij <<
          //         std::endl;

          // truncate PNOs
          size_t pnodrop = 0;
          if (tpno_ != 0.0) {
            for (size_t k = 0; k != occ_ij.rows(); ++k) {
              if (!(occ_ij(k) >= tpno_))
                ++pnodrop;
              else
                break;
            }
          }
          const auto npno = nuocc - pnodrop;

          // add npno to running total of all PNOs in system
          total_pno_sum += npno;

          // Declare eigen matrix pno_trunc
          Eigen::MatrixXd pno_trunc;

          // If npno = 0, substitute a single fake "PNO" with all coefficients
          // equal to zero. All other code will behave as normal
          if (npno == 0) {
            // resize pno_trunc to be size nuocc x 1
            pno_trunc.resize(nuocc_, 1);

            // Create a zero matrix of size nuocc x 1
            Eigen::MatrixXd pno_zero = Eigen::MatrixXd::Zero(nuocc_, 1);

            // Set pno_trunc eqaul to pno_zero matrix
            pno_trunc = pno_zero;
          }  // npno == 0

          // If npno != zero, use actual zet of truncated PNOs
          else {
            pno_trunc.resize(nuocc_, npno);
            pno_trunc = pno_ij.block(0, pnodrop, nuocc, npno);
          }  // npno != 0

          // Store truncated PNOs
          pnos_[i * nocc_act + j] = pno_trunc;

#if PRODUCE_PNO_MOLDEN_FILES
          // write PNOs to Molden
          {
            Eigen::MatrixXd molden_coefs(C_i_eig.rows(), 2 + pno_trunc.cols());
            molden_coefs.col(0) = C_i_eig.col(i);
            molden_coefs.col(1) = C_i_eig.col(j);
            molden_coefs.block(0, 2, C_i_eig.rows(), pno_trunc.cols()) =
                C_a_eig * pno_trunc;

            Eigen::VectorXd occs(2 + pno_trunc.cols());
            occs.setZero();
            occs[0] = 2.0;
            occs[1] = 2.0;

            Eigen::VectorXd evals(2 + pno_trunc.cols());
            evals(0) = 0.0;
            evals(1) = 0.0;
            evals.tail(pno_trunc.cols()) = occ_ij.tail(pno_trunc.cols());

            libint2::molden::Export xport(libint2_atoms, libint2_shells,
                                          molden_coefs, occs, evals);
            xport.write(std::string("pno_") + std::to_string(i) + "_" +
                        std::to_string(j) + ".molden");
          }
#endif

          // Transform F to PNO space
          Eigen::MatrixXd F_pno_ij = pno_trunc.transpose() * F_uocc * pno_trunc;

          // Store just the diagonal elements of F_pno_ij
          // F_pno_diag_[ij] = F_pno_ij.diagonal();
          F_pno_diag_[i * nocc_act + j] = F_pno_ij.diagonal();

          /////// Transform PNOs to canonical PNOs if pno_canonical_ == true

          if (pno_canonical_ && npno > 0) {
            // Compute eigenvectors of F in PNO space
            es.compute(F_pno_ij);
            Eigen::MatrixXd pno_transform_ij = es.eigenvectors();

            // Transform pno_ij to canonical PNO space; pno_ij -> can_pno_ij
            Eigen::MatrixXd can_pno_ij = pno_trunc * pno_transform_ij;

            // Replace standard with canonical PNOs
            // pnos_[ij] = can_pno_ij;
            // F_pno_diag_[ij] = es.eigenvalues();
            pnos_[i * nocc_act + j] = can_pno_ij;
            F_pno_diag_[i * nocc_act + j] = es.eigenvalues();
          }

          // truncate OSVs
          if (i == j) {
            size_t osvdrop = 0;
            if (tosv_ != 0.0) {
              for (size_t k = 0; k != occ_ij.rows(); ++k) {
                if (!(occ_ij(k) >= tosv_))
                  ++osvdrop;
                else
                  break;
              }
            }
            const auto nosv = nuocc - osvdrop;
            if (nosv == 0) {  // all OSV truncated indicates total nonsense
              throw LimitExceeded<size_t>("all OSVs truncated", __FILE__,
                                          __LINE__, 1, 0);
            }

            // Store truncated OSVs
            Eigen::MatrixXd osv_trunc = pno_ij.block(0, osvdrop, nuocc, nosv);
            osvs_[i] = osv_trunc;

            // Transform F to OSV space
            Eigen::MatrixXd F_osv_i =
                osv_trunc.transpose() * F_uocc * osv_trunc;

            // Store just the diagonal elements of F_osv_i
            F_osv_diag_[i] = F_osv_i.diagonal();

            /////// Transform OSVs to canonical OSVs if pno_canonical_ == true
            if (pno_canonical_) {
              // Compute eigenvectors of F in OSV space
              es.compute(F_osv_i);
              Eigen::MatrixXd osv_transform_i = es.eigenvectors();

              // Transform osv_i to canonical OSV space: osv_i -> can_osv_i
              Eigen::MatrixXd can_osv_i = osv_trunc * osv_transform_i;

              // Replace standard with canonical OSVs
              osvs_[i] = can_osv_i;
              F_osv_diag_[i] = es.eigenvalues();
            }
          }
        }
      }
      auto sum_osv = 0;
      for (int i = 0; i < nocc_act; ++i) {
        sum_osv += osvs_[i].cols();
      }
      auto ave_nosv = sum_osv / nocc_act;
      ExEnv::out0() << "The average number of OSVs is " << ave_nosv
                    << std::endl;

      // Compute and print out average number of PNOs per pair
      auto ave_npno = total_pno_sum / (nocc_act * nocc_act);
      ExEnv::out0() << "The average number of PNOs is " << ave_npno
                    << std::endl;

    }  // end if tiling_method_ == rigid

    // When tiling_method_ != rigid

    else {
      auto& fac = factory_;
      auto& world = fac.world();
      auto& ofac = fac.orbital_registry();

      auto nocc = ofac.retrieve("m").rank();
      auto nocc_act = ofac.retrieve("i").rank();
      nocc_act_ = nocc_act;
      auto nuocc = ofac.retrieve("a").rank();
      nuocc_ = nuocc;
      auto nfzc = nocc - nocc_act;

      iter_count_ = 1;

      // Form Fock array
      auto F = fac.compute(L"<p|F|q>[df]");

      // Transform entire Fock array to Eigen Matrix
      Matrix F_all = array_ops::array_to_eigen(F);

      // Select just diagonal elements of Fock aray and transform
      // to Eigen vector; use for computing PNOs
      Vector eps_p = F_all.diagonal();
      auto eps_o = eps_p.segment(nfzc, nocc_act);
      auto eps_v = eps_p.tail(nuocc);

      // Select just the occupied portion of the Fock matrix
      F_occ_act_ = F_all.block(nfzc, nfzc, nocc_act, nocc_act);

      // Select just the unoccupied portion of the Fock matrix
      F_uocc_ = F_all.block(nocc, nocc, nuocc, nuocc);

      // Compute all K_aibj
      // auto K = fac.compute(L"(a b|G|i j)");
      auto K = fac.compute(L"<a b|G|i j>[df]");
      const auto ktrange = K.trange();


      /// Step (1): Convert K to T

      // lambda function to convert K to T; implement using a for_each

      auto form_T = [&eps_o, &eps_v](Tile& result_tile, const Tile& arg_tile) {

        result_tile = Tile(arg_tile.range());

        // determine range of i and j indices
        const int i0 = arg_tile.range().lobound()[2];
        const int in = arg_tile.range().upbound()[2];
        const int j0 = arg_tile.range().lobound()[3];
        const int jn = arg_tile.range().upbound()[3];

        // determine range of a and b indices
        const int a0 = arg_tile.range().lobound()[0];
        const int an = arg_tile.range().upbound()[0];
        const int b0 = arg_tile.range().lobound()[1];
        const int bn = arg_tile.range().upbound()[1];

        auto norm = 0.0;

        // Loop over all four indices to form T
        for (int a = a0, tile_idx = 0; a != an; ++a) {
          const auto e_a = eps_v[a];

          for (int b = b0; b != bn; ++b) {
            const auto e_b = eps_v[b];
            const auto e_ab = e_a + e_b;

            for (int i = i0; i != in; ++i) {
              const auto e_i = eps_o[i];

              for (int j = j0; j != jn; ++j, ++tile_idx) {
                const auto e_j = eps_o[j];
                const auto e_ij = -e_i - e_j;

                const auto e_abij = e_ab + e_ij;
                // const auto e_abij = e_a + e_b - e_i - e_j;
                const auto K_abij = arg_tile[tile_idx];
                const auto T_abij = -K_abij / e_abij;
                const auto abs_result = std::abs(T_abij);
                norm += abs_result * abs_result;
                result_tile[tile_idx] = T_abij;
              }  // j
            }    // i
          }      // b
        }        // a
        return std::sqrt(norm);
      };  // form_T

      auto T2 = TA::foreach(K, form_T);
      world.gop.fence();
      //      std::cout << "Successfully transformed K to T" << std::endl;

      // Reblock T_ so that
      // each occ dim has nocc_act tiles, each of which contains
      // a single element and
      // each uocc dim has single tile containing nuocc elements

      // Create TiledRange1 objects for uocc transformation arrays
      std::vector<std::size_t> uocc_blocks{0, nuocc};

      const TA::TiledRange1 uocc_col =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
      const TA::TiledRange1 uocc_row = ktrange.dim(0);

      // Create TiledRange1 objects for occ transformation arrays
      std::vector<std::size_t> occ_blocks;
      for (std::size_t i = 0; i <= nocc_act_; ++i) {
        occ_blocks.push_back(i);
      }

      const TA::TiledRange1 occ_col =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
      const TA::TiledRange1 occ_row = ktrange.dim(3);

      // Create transition arrays
      reblock_a_ = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      reblock_i_ = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      // Reblock T2
      T T_reblock = detail::reblock_t2(T2, reblock_i_, reblock_a_);

      detail::construct_pno(T_reblock, F_uocc_, tpno_, tosv_,
                            pnos_, npnos_, F_pno_diag_,
                            osvs_, nosvs_, F_osv_diag_, pno_canonical_);

    }  // end if tiling_method_ != rigid

  }

  virtual ~PNOSolver() = default;

  /// @return PNO truncation threshold
  double tpno() const { return tpno_; }
  /// @return OSV truncation threshold
  double tosv() const { return tosv_; }

  const auto& pno(int i, int j) const { return pnos_[i * nocc_act_ + j]; }
  const auto& osv(int i) const { return osvs_[i]; }

  const auto& npnos(int i, int j) const { return npnos_[i * nocc_act_ + j]; }

 private:
  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space

  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {

    T delta_t1_ai;
    T delta_t2_abij;

    // Perform Jacobi update in full space when PNOs being recomputed
    if ((update_pno_ == true) && (iter_count_ % interval_ == 0)) {
      Vector ens_occ_act = F_occ_act_.diagonal();
      Vector ens_uocc = F_uocc_.diagonal();

      delta_t1_ai =
          detail::jacobi_update_t1_ai(r1, ens_occ_act, ens_uocc);
      delta_t2_abij =
          detail::jacobi_update_t2_abij(r2, ens_occ_act, ens_uocc);
    }

    // Perform Jacobi update in PNO space on odd iterations
    else {
      delta_t1_ai =
        detail::pno_jacobi_update_t1(r1, F_occ_act_, F_osv_diag_, osvs_);
      delta_t2_abij =
        detail::pno_jacobi_update_t2(r2, F_occ_act_, F_pno_diag_, pnos_);
    }

    // Reblock delta_t1_ai and delta_t2_abij to match original tiling
    auto delta_t1 = detail::unblock_t1(delta_t1_ai, reblock_i_, reblock_a_);
    auto delta_t2 = detail::unblock_t2(delta_t2_abij, reblock_i_, reblock_a_);
    t1("a,i") += delta_t1("a,i");
    t2("a,b,i,j") += delta_t2("a,b,i,j");
    t1.truncate();
    t2.truncate();

  }

  void update(T& t1, T& t2, const T& r1, const T& r2) override {

    // Reblock r1 and r2
    T r2_reblock = detail::reblock_t2(r2, reblock_i_, reblock_a_);
    T r1_reblock = detail::reblock_t1(r1, reblock_i_, reblock_a_);

    // Whether or not Jacobi update is performed in PNO subspace
    // or full space will be determined within update_only
    update_only(t1, t2, r1_reblock, r2_reblock);

    // Recompute PNOs as appropriate
    if ((update_pno_ == true) && (iter_count_ % interval_ == 0)) {
      T T_reblock = detail::reblock_t2(t2, reblock_i_, reblock_a_);
      detail::reconstruct_pno(T_reblock, F_uocc_, npnos_, nosvs_,
                             pnos_, F_pno_diag_,
                             osvs_, F_osv_diag_, pno_canonical_);
    }



    // transform residuals to the PNO space for the sake of extrapolation
    T r1_osv = detail::osv_transform_ai(r1_reblock, osvs_);
    T r2_pno = detail::pno_transform_abij(r2_reblock, pnos_);
    mpqc::cc::T1T2<T, T> r(r1_osv, r2_pno);
    mpqc::cc::T1T2<T, T> t(t1, t2);
    this->diis().extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
    iter_count_ += 1;
  }

  // squared norm of 1-body residual in OSV subspace
  struct R1SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R1SquaredNormReductionOp(PNOSolver<T, DT>* solver) : solver_(solver) {}

    // Reduction functions
    // Make an empty result object
    result_type operator()() const { return 0; }
    // Post process the result (no operation, passthrough)
    const result_type& operator()(const result_type& result) const {
      return result;
    }
    void operator()(result_type& result, const result_type& arg) const {
      result += arg;
    }
    /// Reduce an argument pair
    void operator()(result_type& result, const argument_type& arg) const {
      const auto i = arg.range().lobound()[1];
      const auto nuocc = arg.range().extent_data()[0];

      const Eigen::MatrixXd arg_osv =
          TA::eigen_map(arg, 1, nuocc) * solver_->osv(i);
      result += arg_osv.squaredNorm();
    }

    PNOSolver<T, DT>* solver_;
  };  // R1SquaredNormReductionOp

  // squared norm of 2-body residual in PNO subspace
  struct R2SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R2SquaredNormReductionOp(PNOSolver<T, DT>* solver) : solver_(solver) {}

    // Reduction functions
    // Make an empty result object
    result_type operator()() const { return 0; }
    // Post process the result (no operation, passthrough)
    const result_type& operator()(const result_type& result) const {
      return result;
    }
    void operator()(result_type& result, const result_type& arg) const {
      result += arg;
    }
    /// Reduce an argument pair
    void operator()(result_type& result, const argument_type& arg) const {
      const auto i = arg.range().lobound()[2];
      const auto j = arg.range().lobound()[3];
      const auto nuocc = arg.range().extent_data()[0];

      const Eigen::MatrixXd arg_pno = solver_->pno(i, j).transpose() *
                                      TA::eigen_map(arg, nuocc, nuocc) *
                                      solver_->pno(i, j);
      result += arg_pno.squaredNorm();
    }

    PNOSolver<T, DT>* solver_;
  };  // R2SquaredNormReductionOp

 public:
  /// Overrides Solver<T,T>::error()
  virtual double error(const T& r1, const T& r2) override {

    auto r1_reblock = detail::reblock_t1(r1, reblock_i_, reblock_a_);
    auto r2_reblock = detail::reblock_t2(r2, reblock_i_, reblock_a_);

    R1SquaredNormReductionOp op1(this);
    R2SquaredNormReductionOp op2(this);

    auto residual = sqrt(r1_reblock("a,i").reduce(op1).get() +
                         r2_reblock("a,b,i,j").reduce(op2).get()) /
                         (size(r1_reblock) + size(r2_reblock));
    if (residual < residual_thresh_) {
      update_pno_ = false;
    }

    return residual;

  }



 private:
  Factory<T, DT>& factory_;
  std::string pno_method_;     //!< the PNO construction method
  bool pno_canonical_;         //!< whether or not to canonicalize PNO/OSV
  bool update_pno_;            //!< whether or not to update PNOs with each CC iter
  std::string tiling_method_;  //!< whether to employ rigid tiling or flexible
                               //! tiling PNO code
  double tpno_;                //!< the PNO truncation threshold
  double tosv_;                //!< the OSV (diagonal PNO) truncation threshold
  int nocc_act_;               //!< the number of active occupied orbitals
  int nuocc_;                  //!< the number of unoccupied orbitals
  T T_;                        //!< the array of MP2 T amplitudes

  int iter_count_;             //!< the CCSD iteration
  int interval_;               //!< the interval at which to update PNOs

  double residual_thresh_;     //!< if residual becomes smaller than this
                               // threshold, stop updating PNOs

  T reblock_i_;
  T reblock_a_;

  Matrix F_occ_act_;
  Matrix F_uocc_;

  // For storing PNOs and and the Fock matrix in the PNO basis
  std::vector<int> npnos_;
  std::vector<Matrix> pnos_;
  std::vector<Vector> F_pno_diag_;


  // For storing OSVs (PNOs when i = j) and the Fock matrix in
  // the OSV basis
  std::vector<int> nosvs_;
  std::vector<Matrix> osvs_;
  std::vector<Vector> F_osv_diag_;



};  // class: PNO solver

////// SVO2 solver ////////

/// SVOSolver updates the CC T amplitudes using standard Jacobi+DIIS in SVO2
/// space
/// @warning This class assumes that the 1- and 2-body amplitudes/residuals
///          given to Solver::update() are laid out as "a,i" and "a,b,i,j",
///          respectively
template <typename T, typename DT>
class SVOSolver : public ::mpqc::cc::DIISSolver<T, T>,
                  public madness::WorldObject<SVOSolver<T, DT>> {
 public:
  typedef typename T::value_type Tile;
  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object; it will be queried for all keywords of ::mpqc::cc::DIISSolver , as well
   * as the following additional keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | pno_method | string | standard | The PNO construction method. Valid values are: \c standard . |
   * | tpno | double | 1e-8 | The PNO construction threshold. This non-negative integer specifies the screening threshold for the eigenvalues of the pair density. Setting this to zero will cause the full (untruncated) set of PNOs to be used. |
   * | tosv | double | 1e-9 | The OSV construction threshold. This non-negative integer specifies the screening threshold for the eigenvalues of the pair density of the diagonal pairs. Setting this to zero will cause the full (untruncated) set of OSVs to be used. |
   */
  // clang-format on
  SVOSolver(const KeyVal& kv, Factory<T, DT>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv),
        madness::WorldObject<SVOSolver<T, DT>>(factory.world()),
        factory_(factory),
        tiling_method_(kv.value<std::string>("tiling_method", "flexible")),
        tsvo2_(kv.value<double>("tsvo2", 1.e-5)),
        tsvo1_(kv.value<double>("tsvo1", 1.e-7)) {
    // part of WorldObject initialization
    this->process_pending();

    // compute and store SVO2s truncated with threshold t_svo2
    // compute and store SVO1s truncated with threshold t_svo1

    // Determine value of tiling_method_
    std::string tiling_method_ = (kv.value<std::string>("tiling_method"));

    // If tiling_method_ = rigid, employ original SVO code
    if (tiling_method_ == "rigid") {
      // Check that tiling is done appropriately
      if (kv.exists("occ_block_size")) {
        int occ_block_size_ = (kv.value<int>("occ_block_size"));
        if (occ_block_size_ != 1) {
          throw InputError(
              "occ_block_size must be set to 1 in the input file.");
        }
      } else {
        throw InputError("occ_block_size was not specified in the input file.");
      }

      if (kv.exists("unocc_block_size")) {
        int unocc_block_size_ = (kv.value<int>("unocc_block_size"));
        if (unocc_block_size_ < 1000000000) {
          throw InputError(
              "unocc_block_size must be greater than or equal to 1000000000 in "
              "the input file.");
        }
      } else {
        throw InputError(
            "unocc_block_size was not specified in the input file.");
      }

      // !!! Original SVO formation code !!! //
      // Do not delete! //

      auto& fac = factory_;
      auto& world = fac.world();
      auto& ofac = fac.orbital_registry();

      auto nocc = ofac.retrieve("m").rank();
      auto nocc_act = ofac.retrieve("i").rank();
      nocc_act_ = nocc_act;
      auto nuocc = ofac.retrieve("a").rank();
      nuocc_ = nuocc;
      auto nfzc = nocc - nocc_act;

      // Form Fock array
      auto F = fac.compute(L"<p|F|q>[df]");

      // Select just diagonal elements of Fock aray and transform
      // to Eigen vector; use for computing SVO2s
      Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
      auto eps_o = eps_p.segment(nfzc, nocc_act);
      auto eps_v = eps_p.tail(nuocc);

      // Transform entire Fock array to Eigen Matrix
      Eigen::MatrixXd F_all = array_ops::array_to_eigen(F);

      // Select just the occupied portion of the Fock matrix
      F_occ_act_ = F_all.block(nfzc, nfzc, nocc_act, nocc_act);

      // Select just the unoccupied portion of the Fock matrix
      Eigen::MatrixXd F_uocc = F_all.block(nocc, nocc, nuocc, nuocc);

      // Compute all K_aibj
      // auto K = fac.compute(L"(a b|G|i j)");
      auto K = fac.compute(L"<a b|G|i j>[df]");
      const auto ktrange = K.trange();

      // zero out amplitudes
      if (!T_.is_initialized()) {
        T_ = T(world, K.trange(), K.shape());
        T_.fill(0.0);
      }

      // For storing SVO2s and and the Fock matrix in the SVO2 basis
      l_svo2s_.resize(nocc_act * nocc_act);
      r_svo2s_.resize(nocc_act * nocc_act);
      F_l_svo2_diag_.resize(nocc_act * nocc_act);
      F_r_svo2_diag_.resize(nocc_act * nocc_act);

      // For storing SVO1s and the Fock matrix in the SVO1 basis
      svo1s_.resize(nocc_act);
      F_svo1_diag_.resize(nocc_act);

      // Keep running sum of all LSVO2s in system for computing average
      auto total_lsvo2_sum = 0;

      // Loop over each pair of occupieds to form amplitude matrices
      for (int i = 0; i < nocc_act; ++i) {
        double eps_i = eps_o[i];

        for (int j = 0; j < nocc_act; ++j) {
          double eps_j = eps_o[j];
          int delta_ij = (i == j) ? 1 : 0;
          std::array<int, 4> tile_ij = {{0, 0, i, j}};
          std::array<int, 4> tile_ji = {{0, 0, j, i}};
          const auto ord_ij = ktrange.tiles_range().ordinal(tile_ij);
          const auto ord_ji = ktrange.tiles_range().ordinal(tile_ji);
          TA::TensorD K_ij = K.find(ord_ij);
          TA::TensorD K_ji = K.find(ord_ji);
          auto ext_ij = K_ij.range().extent_data();
          auto ext_ji = K_ji.range().extent_data();
          Eigen::MatrixXd K_ij_mat =
              TA::eigen_map(K_ij, ext_ij[0] * ext_ij[2], ext_ij[1] * ext_ij[3]);
          Eigen::MatrixXd K_ji_mat =
              TA::eigen_map(K_ji, ext_ji[0] * ext_ji[2], ext_ji[1] * ext_ji[3]);

          Eigen::MatrixXd T_ij(nuocc, nuocc);

          for (int a = 0; a < nuocc; ++a) {
            double eps_a = eps_v[a];
            for (int b = 0; b < nuocc; ++b) {
              double eps_b = eps_v[b];

              T_ij(a, b) = -K_ij_mat(a, b) / (eps_a + eps_b - eps_i - eps_j);

            }  // for each b
          }    // for each a

          // Compute SVD of each T_ij matrix
          Eigen::JacobiSVD<Eigen::MatrixXd> svd(
              T_ij, Eigen::ComputeThinU | Eigen::ComputeThinV);
          auto sing_vals = svd.singularValues();

          // truncate SVO2s
          size_t svo2_keep = 0;
          if (tsvo2_ != 0.0) {
            for (size_t k = 0; k != sing_vals.rows(); ++k) {
              if (sing_vals(k) >= tsvo2_)
                ++svo2_keep;
              else
                break;
            }  // for each k
            //            if (svo2_keep == 0) {  // all SVO2 truncated indicates
            //            total nonsense
            //              throw LimitExceeded<size_t>("all SVO2s truncated",
            //              __FILE__,
            //                                          __LINE__, 1, 0);
            //            }
          }  // if tsvo2 != 0

          else {  // handles case where tsvo2_ == 0
            svo2_keep = nuocc;
          }

          const auto nsvo2 = svo2_keep;
          // std::cout << "nsvo2 = " << nsvo2 << std::endl;

          // add nsvo2 to running total of all LSVO2s in system
          total_lsvo2_sum += nsvo2;

          // Declare eigen matrices to hold truncated SVO2s
          Eigen::MatrixXd r_svo2_trunc;
          Eigen::MatrixXd l_svo2_trunc;

          // If nsvo2 = 0, substitute a single fake "SVO2" with all
          // coefficients equal to zero. All other code will behave
          // as normal

          if (nsvo2 == 0) {
            // resize r_svo2_trunc and l_svo2_trunct to be size nuocc x 1
            r_svo2_trunc.resize(nuocc_, 1);
            l_svo2_trunc.resize(nuocc_, 1);

            // Create a zero matrix of size nuocc x 1
            Eigen::MatrixXd svo2_zero = Eigen::MatrixXd::Zero(nuocc_, 1);

            // Set r_svo2_trunc and l_svo2_trunc equal to svo2_zero matrix
            r_svo2_trunc = svo2_zero;
            l_svo2_trunc = svo2_zero;
          }  // nsvo2 == 0

          // If nsvo2 != 0, use actual set of truncated left and right SVO2s
          else {
            r_svo2_trunc.resize(nuocc_, nsvo2);
            l_svo2_trunc.resize(nuocc_, nsvo2);

            r_svo2_trunc = svd.matrixV().block(0, 0, nuocc, nsvo2);
            l_svo2_trunc = svd.matrixU().block(0, 0, nuocc, nsvo2);
          }  // nsvo2 != 0

          // store truncated SVO2s
          r_svo2s_[i * nocc_act + j] = r_svo2_trunc;
          l_svo2s_[i * nocc_act + j] = l_svo2_trunc;

          // transform F to right SVO2 space and store just diagonal elements
          // Eigen::MatrixXd F_r_svo2 = svd.matrixV().transpose() * F_uocc *
          // svd.matrixV();
          Eigen::MatrixXd F_r_svo2 =
              r_svo2_trunc.transpose() * F_uocc * r_svo2_trunc;
          F_r_svo2_diag_[i * nocc_act + j] = F_r_svo2.diagonal();

          // transform F to left SVO2 space and store just diagonal elements
          // Eigen::MatrixXd F_l_svo2 = svd.matrixU().transpose() *
          // F_uocc.transpose() * svd.matrixU();
          Eigen::MatrixXd F_l_svo2 =
              l_svo2_trunc.transpose() * F_uocc.transpose() * l_svo2_trunc;
          F_l_svo2_diag_[i * nocc_act + j] = F_l_svo2.diagonal();

          // if i == j, truncate and store SVO1s
          if (i == j) {
            size_t svo1_keep = 0;
            if (tsvo1_ != 0.0) {
              for (size_t k = 0; k != sing_vals.rows(); ++k) {
                if (sing_vals(k) >= tsvo1_)
                  ++svo1_keep;
                else
                  break;
              }  // for each k
              if (svo1_keep ==
                  0) {  // all SVO1 truncated indicates total nonsense
                throw LimitExceeded<size_t>("all SVO1s truncated", __FILE__,
                                            __LINE__, 1, 0);
              }
            }  // if tsvo1 != 0

            else {  // handles case where tsvo1_ == 0
              svo1_keep = nuocc;
            }

            const auto nsvo1 = svo1_keep;
            // std::cout << "nsvo1 = " << nsvo1 << std::endl;

            // store truncated SVO1s
            Eigen::MatrixXd svo1_trunc =
                svd.matrixV().block(0, 0, nuocc, nsvo1);
            svo1s_[i] = svo1_trunc;

            // transform F to SVO1 space and store just diagonal elements
            Eigen::MatrixXd F_svo1 =
                svo1_trunc.transpose() * F_uocc * svo1_trunc;
            F_svo1_diag_[i] = F_svo1.diagonal();

          }  // if (i == j)

        }  // for each j
      }    // for each i

      // Compute average number of SVO1s per pair and print out

      auto sum_svo1 = 0;
      for (int i = 0; i < nocc_act; ++i) {
        sum_svo1 += svo1s_[i].cols();
      }
      auto ave_nsvo1 = sum_svo1 / nocc_act;
      // std::cout << "sum_svo1 = " << sum_svo1 << std::endl;
      ExEnv::out0() << "The average number of SVO1s is " << ave_nsvo1
                    << std::endl;

      // Compute average number of SVO2s per pair and print out
      auto ave_nsvo2 = total_lsvo2_sum / (nocc_act_ * nocc_act_);
      ExEnv::out0() << "The average number of SVO2s is " << ave_nsvo2
                    << std::endl;

    }  // if (tiling_method_ == rigid)

    // When tiling_method_ != rigid

    else {
      // Compute integrals

      auto& fac = factory_;
      auto& world = fac.world();
      auto& ofac = fac.orbital_registry();

      auto nocc = ofac.retrieve("m").rank();
      auto nocc_act = ofac.retrieve("i").rank();
      nocc_act_ = nocc_act;
      auto nuocc = ofac.retrieve("a").rank();
      nuocc_ = nuocc;
      auto nfzc = nocc - nocc_act;

      // Form Fock array
      auto F = fac.compute(L"<p|F|q>[df]");

      // Select just diagonal elements of Fock aray and transform
      // to Eigen vector; use for computing SVO2s
      Eigen::VectorXd eps_p = array_ops::array_to_eigen(F).diagonal();
      auto eps_o = eps_p.segment(nfzc, nocc_act);
      auto eps_v = eps_p.tail(nuocc);

      // Transform entire Fock array to Eigen Matrix
      Eigen::MatrixXd F_all = array_ops::array_to_eigen(F);

      // Select just the occupied portion of the Fock matrix
      F_occ_act_ = F_all.block(nfzc, nfzc, nocc_act, nocc_act);

      // Select just the unoccupied portion of the Fock matrix
      Eigen::MatrixXd F_uocc = F_all.block(nocc, nocc, nuocc, nuocc);

      // Compute all K_aibj
      // auto K = fac.compute(L"(a b|G|i j)");
      auto K = fac.compute(L"<a b|G|i j>[df]");
      const auto ktrange = K.trange();

      // Determine number of tiles along each dim of K
      //      const auto ntiles_a = ktrange.dim(0).tile_extent();

      // zero out amplitudes
      if (!T_.is_initialized()) {
        T_ = T(world, K.trange(), K.shape());
        T_.fill(0.0);
      }

      // For storing T_ij matrices
      T_ij_.resize(nocc_act * nocc_act);

      // For storing SVO2s and and the Fock matrix in the SVO2 basis
      l_svo2s_.resize(nocc_act * nocc_act);
      r_svo2s_.resize(nocc_act * nocc_act);
      F_l_svo2_diag_.resize(nocc_act * nocc_act);
      F_r_svo2_diag_.resize(nocc_act * nocc_act);

      // For storing SVO1s and the Fock matrix in the SVO1 basis
      svo1s_.resize(nocc_act);
      F_svo1_diag_.resize(nocc_act);

      /// Step (1): Convert K to T

      // lambda function to convert K to T

      auto form_T = [eps_o, eps_v, this](Tile& result_tile,
                                         const Tile& arg_tile) {

        result_tile = Tile(arg_tile.range());

        // determine range of i and j indices
        const int i0 = arg_tile.range().lobound()[2];
        const int in = arg_tile.range().upbound()[2];
        const int j0 = arg_tile.range().lobound()[3];
        const int jn = arg_tile.range().upbound()[3];

        // determine range of a and b indices
        const int a0 = arg_tile.range().lobound()[0];
        const int an = arg_tile.range().upbound()[0];
        const int b0 = arg_tile.range().lobound()[1];
        const int bn = arg_tile.range().upbound()[1];

        auto norm = 0.0;

        // Loop over all four indices to form T
        for (int a = a0, tile_idx = 0; a != an; ++a) {
          const auto e_a = eps_v[a];

          for (int b = b0; b != bn; ++b) {
            const auto e_b = eps_v[b];
            const auto e_ab = e_a + e_b;

            for (int i = i0; i != in; ++i) {
              const auto e_i = eps_o[i];

              for (int j = j0; j != jn; ++j, ++tile_idx) {
                const auto e_j = eps_o[j];
                const auto e_ij = -e_i - e_j;

                const auto e_abij = e_ab + e_ij;
                // const auto e_abij = e_a + e_b - e_i - e_j;
                const auto K_abij = arg_tile[tile_idx];
                const auto T_abij = -K_abij / e_abij;
                const auto abs_result = std::abs(T_abij);
                norm += abs_result * abs_result;
                result_tile[tile_idx] = T_abij;
              }  // j
            }    // i
          }      // b
        }        // a
        return std::sqrt(norm);
      };  // form_T

      auto T_ = TA::foreach (K, form_T);
      T_.world().gop.fence();
      std::cout << "Successfully transformed K to T" << std::endl;

      // Reblock T_ so that
      // each occ dim has nocc_act tiles, each of which contains
      // a single element and
      // each uocc dim has single tile containing nuocc elements

      // Create TiledRange1 objects for uocc transformation arrays
      std::vector<std::size_t> uocc_blocks{0, nuocc};

      const TA::TiledRange1 uocc_col =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
      const TA::TiledRange1 uocc_row = ktrange.dim(0);

      // Create TiledRange1 objects for occ transformation arrays
      std::vector<std::size_t> occ_blocks;
      for (std::size_t i = 0; i <= nocc_act; ++i) {
        occ_blocks.push_back(i);
      }

      const TA::TiledRange1 occ_col =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
      const TA::TiledRange1 occ_row = ktrange.dim(3);

      std::cout << "uocc_row number of tiles: " << uocc_row.tile_extent()
                << std::endl;
      std::cout << "uocc_col number of tiles: " << uocc_col.tile_extent()
                << std::endl;

      std::cout << "uocc_row number of elements: " << uocc_row.extent()
                << std::endl;
      std::cout << "uocc_col number of elements: " << uocc_col.extent()
                << std::endl;

      std::cout << "occ_row number of tiles: " << occ_row.tile_extent()
                << std::endl;
      std::cout << "occ_col number of tiles: " << occ_col.tile_extent()
                << std::endl;

      std::cout << "occ_row number of elements: " << occ_row.extent()
                << std::endl;
      std::cout << "occ_col number of elements: " << occ_col.extent()
                << std::endl;

      // Create transition arrays
      T a_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      T b_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      T i_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      T j_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      // Reblock T_
      T T_reblock;

      T_reblock("an,bn,in,jn") = T_("a,b,i,j") * j_trans_array("j,jn") *
                                 i_trans_array("i,in") * b_trans_array("b,bn") *
                                 a_trans_array("a,an");

      std::cout << "reblocking step worked" << std::endl;
      std::cout << "T_trange:\n" << T_.trange() << std::endl;
      std::cout << "Reblocked T_ trange:\n" << T_reblock.trange() << std::endl;

      // Transform each tile of T_reblock to an eigen matrix and store in T_ij_
      for (auto i = 0; i != nocc_act_; ++i) {
        for (auto j = 0; j != nocc_act_; ++j) {
          auto ij = i * nocc_act_ + j;
          Tile T_ij_tile = T_reblock.find({0, 0, i, j});
          T_ij_[ij] = TA::eigen_map(T_ij_tile, nuocc_, nuocc_);

        }  // j
      }    // i

      // Keep running sum of all LSVO2s in system for computing average
      auto total_lsvo2_sum = 0;

      // For each i,j pair, compute SVD of T_ij, form and truncate SVOs, store
      // SVOs
      for (auto i = 0; i != nocc_act_; ++i) {
        for (auto j = 0; j != nocc_act_; ++j) {
          auto ij = i * nocc_act_ + j;
          Eigen::MatrixXd T_ij = T_ij_[ij];

          // Compute SVD of each T_ij matrix
          Eigen::JacobiSVD<Eigen::MatrixXd> svd(
              T_ij, Eigen::ComputeThinU | Eigen::ComputeThinV);
          auto sing_vals = svd.singularValues();

          // truncate SVO2s
          size_t svo2_keep = 0;
          if (tsvo2_ != 0.0) {
            for (size_t k = 0; k != sing_vals.rows(); ++k) {
              if (sing_vals(k) >= tsvo2_)
                ++svo2_keep;
              else
                break;
            }  // for each k

            //            if (svo2_keep == 0) {  // all SVO2 truncated indicates
            //            total nonsense
            //              throw LimitExceeded<size_t>("all SVO2s truncated",
            //              __FILE__,
            //                                          __LINE__, 1, 0);
            //            }
          }  // if tsvo2 != 0

          else {  // handles case where tsvo2_ == 0
            svo2_keep = nuocc;
          }

          const auto nsvo2 = svo2_keep;

          // add nosvo2 to running total of all LSVO2s in system
          total_lsvo2_sum += nsvo2;

          // Declare eigen matrices to hold truncated SVO2s
          Eigen::MatrixXd r_svo2_trunc;
          Eigen::MatrixXd l_svo2_trunc;

          // If nsvo2 = 0. substitute a single fake "SVO2" with all
          // coefficients equal to zero. All other code will behave
          // as normal

          if (nsvo2 == 0) {
            // resize r_svo2_trunc and l_svo2_trunct to be size nuocc x 1
            r_svo2_trunc.resize(nuocc_, 1);
            l_svo2_trunc.resize(nuocc_, 1);

            // Create a zero matrix of size nuocc x 1
            Eigen::MatrixXd svo2_zero = Eigen::MatrixXd::Zero(nuocc_, 1);

            // Set r_svo2_trunc and l_svo2_trunc equal to svo2_zero matrix
            r_svo2_trunc = svo2_zero;
            l_svo2_trunc = svo2_zero;
          }  // nsvo2 == 0

          // If nsvo2 != 0, use actual set of truncated left and right SVO2s
          else {
            r_svo2_trunc.resize(nuocc_, nsvo2);
            l_svo2_trunc.resize(nuocc_, nsvo2);

            r_svo2_trunc = svd.matrixV().block(0, 0, nuocc, nsvo2);
            l_svo2_trunc = svd.matrixU().block(0, 0, nuocc, nsvo2);
          }  // nsvo2 != 0

          // store truncated SVO2s
          r_svo2s_[i * nocc_act + j] = r_svo2_trunc;
          l_svo2s_[i * nocc_act + j] = l_svo2_trunc;

          // transform F to right SVO2 space and store just diagonal elements
          Eigen::MatrixXd F_r_svo2 =
              r_svo2_trunc.transpose() * F_uocc * r_svo2_trunc;
          F_r_svo2_diag_[ij] = F_r_svo2.diagonal();

          // transform F to left SVO2 space and store just diagonal elements
          Eigen::MatrixXd F_l_svo2 =
              l_svo2_trunc.transpose() * F_uocc.transpose() * l_svo2_trunc;
          F_l_svo2_diag_[ij] = F_l_svo2.diagonal();

          // if i == j, truncate and store SVO1s
          if (i == j) {
            size_t svo1_keep = 0;
            if (tsvo1_ != 0.0) {
              for (size_t k = 0; k != sing_vals.rows(); ++k) {
                if (sing_vals(k) >= tsvo1_)
                  ++svo1_keep;
                else
                  break;
              }  // for each k
              if (svo1_keep ==
                  0) {  // all SVO1 truncated indicates total nonsense
                throw LimitExceeded<size_t>("all SVO1s truncated", __FILE__,
                                            __LINE__, 1, 0);
              }
            }  // if tsvo1 != 0

            else {  // handles case where tsvo1_ == 0
              svo1_keep = nuocc;
            }

            const auto nsvo1 = svo1_keep;
            // std::cout << "nsvo1 = " << nsvo1 << std::endl;

            // store truncated SVO1s
            Eigen::MatrixXd svo1_trunc =
                svd.matrixV().block(0, 0, nuocc, nsvo1);
            svo1s_[i] = svo1_trunc;

            // transform F to SVO1 space and store just diagonal elements
            Eigen::MatrixXd F_svo1 =
                svo1_trunc.transpose() * F_uocc * svo1_trunc;
            F_svo1_diag_[i] = F_svo1.diagonal();

          }  // if (i == j)

        }  // for each j
      }    // for each i

      // Compute average number of SVO1s per pair and print out

      auto sum_svo1 = 0;
      for (int i = 0; i < nocc_act; ++i) {
        sum_svo1 += svo1s_[i].cols();
      }
      auto ave_nsvo1 = sum_svo1 / nocc_act;
      // std::cout << "sum_svo1 = " << sum_svo1 << std::endl;
      ExEnv::out0() << "The average number of SVO1s is " << ave_nsvo1
                    << std::endl;

      // Compute average number of SVO2s per pair and print out
      auto ave_nsvo2 = total_lsvo2_sum / (nocc_act_ * nocc_act_);
      ExEnv::out0() << "The average number of SVO2s is " << ave_nsvo2
                    << std::endl;

    }  // tilinig_method_ != rigid

  }  // SVOSolver

  virtual ~SVOSolver() = default;

  /// @return SVO2 truncation threshold
  double tsvo2() const { return tsvo2_; }
  /// @return SVO1 truncation threshold
  double tsvo1() const { return tsvo1_; }

  const auto& l_svo2(int i, int j) const { return l_svo2s_[i * nocc_act_ + j]; }
  const auto& r_svo2(int i, int j) const { return r_svo2s_[i * nocc_act_ + j]; }
  const auto& sov1(int i) const { return svo1s_[i]; }

 private:
  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space
  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    auto delta_t1_ai = jacobi_update_t1(r1, F_occ_act_, F_svo1_diag_, svo1s_);
    auto delta_t2_abij = jacobi_update_t2(r2, F_occ_act_, F_l_svo2_diag_,
                                          F_r_svo2_diag_, l_svo2s_, r_svo2s_);

    // Reblock delta_t1_ai and delta_t2_abij to match original tiling
    auto delta_t1 = unblock_delta_t1(delta_t1_ai);
    auto delta_t2 = unblock_delta_t2(delta_t2_abij);
    t1("a,i") += delta_t1("a,i");
    t2("a,b,i,j") += delta_t2("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }

  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    // reblock r1 and r2
    T r2_reblock = reblock_r2(r2);
    T r1_reblock = reblock_r1(r1);

    update_only(t1, t2, r1_reblock, r2_reblock);
    T r1_svo1 = svo_transform_ai(r1_reblock, svo1s_);
    T r2_svo2 = svo_transform_abij(r2_reblock, l_svo2s_, r_svo2s_);
    mpqc::cc::T1T2<T, T> r(r1_svo1, r2_svo2);
    mpqc::cc::T1T2<T, T> t(t1, t2);
    this->diis().extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> reblock_r2(
      const TA::DistArray<Tile, Policy>& r2) {
    auto reblock2 = [this](const T& arg) {
      auto T_trange = T_.trange();
      auto& world = factory_.world();

      // Create TiledRange1 objects for uocc transformation arrays
      std::vector<int> uocc_blocks{0, nuocc_};

      const TA::TiledRange1 uocc_col =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
      const TA::TiledRange1 uocc_row = T_trange.dim(0);

      // Create TiledRange1 objects for occ transformation arrays
      std::vector<int> occ_blocks;
      for (std::size_t i = 0; i <= nocc_act_; ++i) {
        occ_blocks.push_back(i);
      }

      const TA::TiledRange1 occ_col =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
      const TA::TiledRange1 occ_row = T_trange.dim(3);

      // Create occ and uocc transformation arrays
      T uocc_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      T occ_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      // Reblock r2
      T result;
      result("an,bn,in,jn") =
          arg("a,b,i,j") * occ_trans_array("j,jn") * occ_trans_array("i,in") *
          uocc_trans_array("b,bn") * uocc_trans_array("a,an");

      return result;

    };
    auto r2_reblock = reblock2(r2);
    return r2_reblock;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> reblock_r1(
      const TA::DistArray<Tile, Policy>& r1) {
    auto reblock1 = [this](const T& arg) {
      auto T_trange = T_.trange();
      auto& world = factory_.world();

      // Create TiledRange1 objects for uocc transformation arrays
      std::vector<int> uocc_blocks{0, nuocc_};

      const TA::TiledRange1 uocc_col =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
      const TA::TiledRange1 uocc_row = T_trange.dim(0);

      // Create TiledRange1 objects for occ transformation arrays
      std::vector<int> occ_blocks;
      for (std::size_t i = 0; i <= nocc_act_; ++i) {
        occ_blocks.push_back(i);
      }

      const TA::TiledRange1 occ_col =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
      const TA::TiledRange1 occ_row = T_trange.dim(3);

      // Create occ and uocc transformation arrays
      T uocc_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      T occ_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      // Reblock r1
      T result;
      result("an,in") =
          arg("a,i") * occ_trans_array("i,in") * uocc_trans_array("a,an");

      return result;

    };
    auto r1_reblock = reblock1(r1);
    return r1_reblock;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> unblock_delta_t2(
      const TA::DistArray<Tile, Policy>& delta_t2_abij) {
    auto unblock2 = [this](const T& arg) {
      auto T_trange = T_.trange();
      auto& world = factory_.world();

      // Create TiledRange1 objects for uocc transformation arrays
      std::vector<int> uocc_blocks{0, nuocc_};

      const TA::TiledRange1 uocc_row =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
      const TA::TiledRange1 uocc_col = T_trange.dim(0);

      // Create TiledRange1 objects for occ transformation arrays
      std::vector<int> occ_blocks;
      for (std::size_t i = 0; i <= nocc_act_; ++i) {
        occ_blocks.push_back(i);
      }

      const TA::TiledRange1 occ_row =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
      const TA::TiledRange1 occ_col = T_trange.dim(3);

      // Create occ and uocc transformation arrays
      T uocc_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      T occ_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      // Reblock delta_t2_abij to original blocking
      T result;
      result("an,bn,in,jn") =
          arg("a,b,i,j") * occ_trans_array("j,jn") * occ_trans_array("i,in") *
          uocc_trans_array("b,bn") * uocc_trans_array("a,an");

      return result;

    };
    auto delta_t2 = unblock2(delta_t2_abij);
    return delta_t2;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> unblock_delta_t1(
      const TA::DistArray<Tile, Policy>& delta_t1_ai) {
    auto unblock1 = [this](const T& arg) {
      auto T_trange = T_.trange();
      auto& world = factory_.world();

      // Create TiledRange1 objects for uocc transformation arrays
      std::vector<int> uocc_blocks{0, nuocc_};

      const TA::TiledRange1 uocc_row =
          TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
      const TA::TiledRange1 uocc_col = T_trange.dim(0);

      // Create TiledRange1 objects for occ transformation arrays
      std::vector<int> occ_blocks;
      for (std::size_t i = 0; i <= nocc_act_; ++i) {
        occ_blocks.push_back(i);
      }

      const TA::TiledRange1 occ_row =
          TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
      const TA::TiledRange1 occ_col = T_trange.dim(3);

      // Create occ and uocc transformation arrays
      T uocc_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

      T occ_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<
          Tile, TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

      // Reblock delta_t1_ai to original blocking
      T result;
      result("an,in") =
          arg("a,i") * occ_trans_array("i,in") * uocc_trans_array("a,an");

      return result;

    };
    auto delta_t1 = unblock1(delta_t1_ai);
    return delta_t1;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t2(
      const TA::DistArray<Tile, Policy>& r2_abij,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_l_svo2_diag,
      const std::vector<Eigen::VectorXd>& F_r_svo2_diag,
      const std::vector<Eigen::MatrixXd>& l_svo2s,
      const std::vector<Eigen::MatrixXd>& r_svo2s) {
    auto update2 = [F_occ_act, F_l_svo2_diag, F_r_svo2_diag, l_svo2s, r_svo2s,
                    this](Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd l_svo2_ij = l_svo2s[ij];
      Eigen::MatrixXd r_svo2_ij = r_svo2s[ij];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to SVO2 basis
      const Eigen::MatrixXd r2_svo2 =
          l_svo2_ij.transpose() *
          TA::eigen_map(arg_tile, ext[0] * ext[2], ext[1] * ext[3]) * r_svo2_ij;

      // Create a matrix delta_t2_pno to hold updated values of delta_t2 in PNO
      // basis this matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::MatrixXd delta_t2_svo2 = r2_svo2;

      // Select correct vector containing diagonal elements of Fock matrix in
      // PNO basis
      const Eigen::VectorXd& l_uocc = F_l_svo2_diag[ij];
      const Eigen::VectorXd& r_uocc = F_r_svo2_diag[ij];

      // Determine number of PNOs
      const auto npno = l_uocc.rows();

      // Determine number of uocc
      const auto nuocc = l_svo2_ij.rows();

      // Select e_i and e_j
      const auto e_i = F_occ_act(i, i);
      const auto e_j = F_occ_act(j, j);

      for (auto a = 0; a < npno; ++a) {
        const auto e_a = l_uocc[a];
        for (auto b = 0; b < npno; ++b) {
          const auto e_b = r_uocc[b];
          const auto e_abij = e_i + e_j - e_a - e_b;
          const auto r_abij = r2_svo2(a, b);
          delta_t2_svo2(a, b) = r_abij / e_abij;
        }
      }

      // Back transform delta_t2_svo2 to full space
      // Eigen::MatrixXd delta_t2_full =
      //     svo1_ij * delta_t2_svo2 * l_svo2_ij.transpose();

      Eigen::MatrixXd delta_t2_full =
          l_svo2_ij * delta_t2_svo2 * r_svo2_ij.transpose();

      // Convert delta_t2_full to tile and compute norm
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nuocc; ++r) {
        for (auto c = 0; c < nuocc; ++c) {
          const auto idx = r * nuocc + c;
          const auto elem = delta_t2_full(r, c);
          const auto abs_elem = std::abs(elem);
          norm += abs_elem * abs_elem;
          result_tile[idx] = elem;
        }
      }

      return std::sqrt(norm);
    };

    auto delta_t2_abij = TA::foreach (r2_abij, update2);
    delta_t2_abij.world().gop.fence();
    return delta_t2_abij;
  }  // jacobi_update_t2

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t1(
      const TA::DistArray<Tile, Policy>& r1_ai,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_svo1_diag,
      const std::vector<Eigen::MatrixXd>& svo1s) {
    auto update1 = [F_occ_act, F_svo1_diag, svo1s, this](Tile& result_tile,
                                                         const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i index
      const auto i = arg_tile.range().lobound()[1];

      // Select appropriate matrix of SVO1s
      Eigen::MatrixXd svo1_i = svo1s[i];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to SVO2 basis
      const Eigen::VectorXd r1_svo1 =
          svo1_i.transpose() * TA::eigen_map(arg_tile, ext[0], ext[1]);

      // Create a matrix delta_t1_svo1 to hold updated values of delta t1 in SVO
      // basis. This matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::VectorXd delta_t1_svo1 = r1_svo1;

      // Select correct vector containing diagonal elements of Fock matrix in
      // SVO2 basis
      const Eigen::VectorXd& uocc = F_svo1_diag[i];

      // Determine number of SVO1s
      const auto nsvo1 = uocc.rows();

      // Determine number of uocc
      const auto nuocc = svo1_i.rows();

      // Select e_i
      const auto e_i = F_occ_act(i, i);

      for (auto a = 0; a < nsvo1; ++a) {
        const auto e_a = uocc[a];
        const auto e_ai = e_i - e_a;
        const auto r_ai = r1_svo1(a);
        delta_t1_svo1(a) = r_ai / e_ai;
      }

      // Back transform delta_t1_svo1 to full space
      // Eigen::MatrixXd delta_t1_full = osv_i * delta_t1_osv *
      // osv_i.transpose();
      Eigen::VectorXd delta_t1_full = svo1_i * delta_t1_svo1;
      // Eigen::VectorXd delta_t1_full = delta_t1_svo1 * svo1_i;

      // Convert delta_t1_full to tile and compute norm
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nuocc; ++r) {
        const auto elem = delta_t1_full(r);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto delta_t1_ai = TA::foreach (r1_ai, update1);
    delta_t1_ai.world().gop.fence();
    return delta_t1_ai;
  }  // jacobi_update_t1

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> svo_transform_abij(
      const TA::DistArray<Tile, Policy>& abij,
      const std::vector<Eigen::MatrixXd>& l_svo2s,
      const std::vector<Eigen::MatrixXd>& r_svo2s) {
    auto tform = [l_svo2s, r_svo2s, this](Tile& result_tile,
                                          const Tile& arg_tile) {

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of SVO2s
      const auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd l_svo2_ij = l_svo2s[ij];
      Eigen::MatrixXd r_svo2_ij = r_svo2s[ij];
      const auto nuocc = l_svo2_ij.rows();
      const auto nsvo2 = l_svo2_ij.cols();

      // Convert data in tile to Eigen::Map and transform to SVO2 basis
      const Eigen::MatrixXd result_eig = l_svo2_ij.transpose() *
                                         TA::eigen_map(arg_tile, nuocc, nuocc) *
                                         r_svo2_ij;

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{nsvo2, nsvo2, 1l, 1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nsvo2; ++r) {
        for (auto c = 0; c < nsvo2; ++c) {
          const auto idx = r * nsvo2 + c;
          const auto elem = result_eig(r, c);
          const auto abs_elem = std::abs(elem);
          norm += abs_elem * abs_elem;
          result_tile[idx] = elem;
        }
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach (abij, tform);
    result.world().gop.fence();
    return result;
  }  // svo_transform_abij

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> svo_transform_ai(
      const TA::DistArray<Tile, Policy>& ai,
      const std::vector<Eigen::MatrixXd>& svo1s) {
    auto tform = [svo1s, this](Tile& result_tile, const Tile& arg_tile) {

      // determine i index
      const auto i = arg_tile.range().lobound()[1];

      // Select appropriate matrix of SVOs
      Eigen::MatrixXd svo1_i = svo1s[i];
      const auto nuocc = svo1_i.rows();
      const auto nsvo1 = svo1_i.cols();

      // Convert data in tile to Eigen::Map and transform to SVO1 basis
      const Eigen::MatrixXd result_eig =
          svo1_i.transpose() * TA::eigen_map(arg_tile, nuocc, 1);

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{nsvo1, 1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nsvo1; ++r) {
        const auto elem = result_eig(r, 0);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach (ai, tform);
    result.world().gop.fence();
    return result;
  }  // svo_transform_ai

  // squared norm of 1-body residual in SVO1 subspace
  struct R1SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R1SquaredNormReductionOp(SVOSolver<T, DT>* solver) : solver_(solver) {}

    // Reduction functions
    // Make an empty result object
    result_type operator()() const { return 0; }
    // Post process the result (no operation, passthrough)
    const result_type& operator()(const result_type& result) const {
      return result;
    }
    void operator()(result_type& result, const result_type& arg) const {
      result += arg;
    }
    /// Reduce an argument pair
    void operator()(result_type& result, const argument_type& arg) const {
      const auto i = arg.range().lobound()[1];
      const auto nuocc = arg.range().extent_data()[0];
      const Eigen::MatrixXd arg_svo1 =
          TA::eigen_map(arg, 1, nuocc) * solver_->svo1s_[i];
      result += arg_svo1.squaredNorm();
    }

    SVOSolver<T, DT>* solver_;
  };  // R1SquaredNormReductionOp

  // squared norm of 2-body residual in SVO2 subspace
  struct R2SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R2SquaredNormReductionOp(SVOSolver<T, DT>* solver) : solver_(solver) {}

    // Reduction functions
    // Make an empty result object
    result_type operator()() const { return 0; }
    // Post process the result (no operation, passthrough)
    const result_type& operator()(const result_type& result) const {
      return result;
    }
    void operator()(result_type& result, const result_type& arg) const {
      result += arg;
    }
    /// Reduce an argument pair
    void operator()(result_type& result, const argument_type& arg) const {
      const auto i = arg.range().lobound()[2];
      const auto j = arg.range().lobound()[3];
      const auto nuocc = arg.range().extent_data()[0];
      const Eigen::MatrixXd arg_svo2 = solver_->l_svo2(i, j).transpose() *
                                       TA::eigen_map(arg, nuocc, nuocc) *
                                       solver_->r_svo2(i, j);
      result += arg_svo2.squaredNorm();
    }

    SVOSolver<T, DT>* solver_;
  };  // R2SquaredNormReductionOp

 public:
  /// Overrides Solver<T,T>::error()
  virtual double error(const T& r1, const T& r2) override {
    auto r1_reblock = reblock_r1(r1);
    auto r2_reblock = reblock_r2(r2);

    R1SquaredNormReductionOp op1(this);
    R2SquaredNormReductionOp op2(this);
    return sqrt(r1_reblock("a,i").reduce(op1).get() +
                r2_reblock("a,b,i,j").reduce(op2).get()) /
           (size(r1_reblock) + size(r2_reblock));
  }

 private:
  Factory<T, DT>& factory_;
  std::string tiling_method_;  //!< whether to employ rigid tiling or flexible
                               //!< tiling SVO solver
  double tsvo2_;               //!< the truncation threshold for SVO2s
  double tsvo1_;               //!< the truncation threshold for SVO1s
  int nocc_act_;               //!< the number of active occupied orbitals
  int nuocc_;                  //!< the number of unoccupied orbitals
  T T_;

  Eigen::MatrixXd F_occ_act_;

  // For storing T_ij matrices
  std::vector<Eigen::MatrixXd> T_ij_;

  // For storing SVO2s and the Fock matrix in the SVO2 basis
  std::vector<Eigen::MatrixXd> l_svo2s_;
  std::vector<Eigen::MatrixXd> r_svo2s_;
  std::vector<Eigen::VectorXd> F_l_svo2_diag_;
  std::vector<Eigen::VectorXd> F_r_svo2_diag_;

  // For storing SVO1s and the Fock matrix in the SVO1 basis
  std::vector<Eigen::MatrixXd> svo1s_;
  std::vector<Eigen::VectorXd> F_svo1_diag_;

};  // class: SVOSolver

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_ */
