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
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/math/linalg/diagonal_array.h"

using Array = TA::DistArray<TA::TensorD, TA::SparsePolicy>;

namespace mpqc {
namespace lcao {
namespace cc {

namespace detail {

template <typename Tile, typename Policy,
          typename EigenVectorX =
              Eigen::Matrix<typename Tile::element_type, Eigen::Dynamic, 1>>
TA::DistArray<Tile, Policy> jacobi_update_t2_abij(
    const TA::DistArray<Tile, Policy>& r2_abij, const EigenVectorX& ens_occ,
    const EigenVectorX& ens_uocc) {
  auto denom = [ens_occ, ens_uocc](Tile& result_tile, const Tile& arg_tile) {

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

template <typename Tile, typename Policy,
          typename EigenVectorX =
              Eigen::Matrix<typename Tile::element_type, Eigen::Dynamic, 1>>
TA::DistArray<Tile, Policy> jacobi_update_t1_ai(
    const TA::DistArray<Tile, Policy>& r1_ai, const EigenVectorX& ens_occ,
    const EigenVectorX& ens_uocc) {
  auto denom = [ens_occ, ens_uocc](Tile& result_tile, const Tile& arg_tile) {

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
template <typename T>
// template <typename Tile, typename Policy>
class PNOSolver : public ::mpqc::cc::DIISSolver<T, T>,
                  public madness::WorldObject<PNOSolver<T>> {
 public:

  typedef typename T::value_type Tile;
//  typedef typename T::

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
  PNOSolver(const KeyVal& kv, Factory<T>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv),
        madness::WorldObject<PNOSolver<T>>(factory.world()),
        factory_(factory),
        pno_method_(kv.value<std::string>("pno_method", "standard")),
        pno_canonical_(kv.value<std::string>("pno_canonical", "false")),
        tpno_(kv.value<double>("tpno", 1.e-8)),
        tosv_(kv.value<double>("tosv", 1.e-9)) {
    // part of WorldObject initialization
    this->process_pending();

    // compute and store PNOs truncated with threshold tpno_
    // store PNOs for diagonal pair as OSVs truncated with threshold tosv_

    // // Check that tiling is done appropriately
    // if (kv.exists("occ_block_size")) {
    //   int occ_block_size_ = (kv.value<int>("occ_block_size"));
    //   if (occ_block_size_ != 1) {
    //     throw InputError("occ_block_size must be set to 1 in the input file.");
    //   }
    // } else {
    //   throw InputError("occ_block_size was not specified in the input file.");
    // }

    // if (kv.exists("unocc_block_size")) {
    //   int unocc_block_size_ = (kv.value<int>("unocc_block_size"));
    //   if (unocc_block_size_ < 1000000000) {
    //     throw InputError(
    //         "unocc_block_size must be greater than or equal to 1000000000 in "
    //         "the input file.");
    //   }
    // } else {
    //   throw InputError("unocc_block_size was not specified in the input file.");
    // }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

    auto& fac = factory_;
    auto& world = fac.world();
    auto& ofac = fac.orbital_registry();

    auto nocc = ofac.retrieve("m").rank();
    auto nocc_act = ofac.retrieve("i").rank();
    nocc_act_ = nocc_act;
    auto nuocc = ofac.retrieve("a").rank();
    auto nfzc = nocc - nocc_act;

    // Form Fock array
    auto F = fac.compute(L"<p|F|q>[df]");

    // Select just diagonal elements of Fock aray and transform
    // to Eigen vector; use for computing PNOs
    Eigen::VectorXd eps_p = TA::array_to_eigen(F).diagonal();
    auto eps_o = eps_p.segment(nfzc, nocc_act);
    auto eps_v = eps_p.tail(nuocc);

    // Transform entire Fock array to Eigen Matrix
    Eigen::MatrixXd F_all = TA::array_to_eigen(F);

    // Select just the occupied portion of the Fock matrix
    F_occ_act_ = F_all.block(nfzc, nfzc, nocc_act, nocc_act);

    // Select just the unoccupied portion of the Fock matrix
    Eigen::MatrixXd F_uocc = F_all.block(nocc, nocc, nuocc, nuocc);

    // Compute all K_aibj
    // auto K = fac.compute(L"(a b|G|i j)");
    auto K = fac.compute(L"<a b|G|i j>[df]");
    const auto ktrange = K.trange();

    // Determine number of tiles along each dim of K
    const auto ntiles_a = ktrange.dim(0).tile_extent();
    const auto ntiles_b = ktrange.dim(1).tile_extent();
    const auto ntiles_i = ktrange.dim(2).tile_extent();
    const auto ntiles_j = ktrange.dim(3).tile_extent();

    // zero out amplitudes
    if (!T_.is_initialized()) {
      T_ = Array(world, K.trange(), K.shape());
      T_.fill(0.0);
    }

    // For storing PNOs and and the Fock matrix in the PNO basis
    pnos_.resize(nocc_act * nocc_act);
    F_pno_diag_.resize(nocc_act * nocc_act);

    // For storing OSVs (PNOs when i = j) and the Fock matrix in
    // the OSV basis
    osvs_.resize(nocc_act);
    F_osv_diag_.resize(nocc_act);




    /// Step (1): Convert K to T

    // lambda function to convert K to T

    auto form_T = [eps_o, eps_v, this](
                     Tile& result_tile, const Tile& arg_tile) {


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

        for (int b = b0; b!= bn; ++b) {
          const auto e_b = eps_v[b];
          const auto e_ab = e_a + e_b;

          for (int i = i0; i != in; ++i) {
            const auto e_i = eps_o[i];

            for (int j = j0; j != jn; ++j, ++tile_idx) {
              const auto e_j = eps_o[j];
              const auto e_ij = -e_i - e_j;

              const auto e_abij = e_ab + e_ij;
              //const auto e_abij = e_a + e_b - e_i - e_j;
              const auto K_abij = arg_tile[tile_idx];
              const auto T_abij = -K_abij / e_abij;
              const auto abs_result = std::abs(T_abij);
              norm += abs_result * abs_result;
              result_tile[tile_idx] = T_abij;
            } // j
          } // i
        } // b
      } // a
      return std::sqrt(norm);
    };

    auto T_ = TA::foreach(K, form_T);
    T_.world().gop.fence();
    std::cout << "Successfully transformed K to T" << std::endl;


    // Reblock T_ so that
    // each occ dim has nocc_act tiles, each of which contains
    // a single element and
    // each uocc dim has single tile containing nuocc elements


    // Create TiledRange1 objects for uocc transformation arrays
    std::vector<std::size_t> uocc_blocks {0, nuocc};

    const TA::TiledRange1 uocc_col = TA::TiledRange1(uocc_blocks.begin(), uocc_blocks.end());
    const TA::TiledRange1 uocc_row = ktrange.dim(0);


    // Create TiledRange1 objects for occ transformation arrays
    std::vector<std::size_t> occ_blocks;
    for (std::size_t i = 0; i <= nocc_act; ++i) {
        occ_blocks.push_back(i);
    }

    const TA::TiledRange1 occ_col = TA::TiledRange1(occ_blocks.begin(), occ_blocks.end());
    const TA::TiledRange1 occ_row = ktrange.dim(3);

    std::cout << "uocc_row number of tiles: " << uocc_row.tile_extent() << std::endl;
    std::cout << "uocc_col number of tiles: " << uocc_col.tile_extent() << std::endl;

    std::cout << "uocc_row number of elements: " << uocc_row.extent() << std::endl;
    std::cout << "uocc_col number of elements: " << uocc_col.extent() << std::endl;

    std::cout << "occ_row number of tiles: " << occ_row.tile_extent() << std::endl;
    std::cout << "occ_col number of tiles: " << occ_col.tile_extent() << std::endl;

    std::cout << "occ_row number of elements: " << occ_row.extent() << std::endl;
    std::cout << "occ_col number of elements: " << occ_col.extent() << std::endl;


    // Create transition arrays
    T a_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<Tile,
                            TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

    T b_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<Tile,
                            TA::detail::policy_t<T>>(world, uocc_row, uocc_col, 1.0);

    T i_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<Tile,
                            TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);

    T j_trans_array = mpqc::array_ops::create_diagonal_array_from_eigen<Tile,
                            TA::detail::policy_t<T>>(world, occ_row, occ_col, 1.0);



    // Reblock T_
    T T_reblock;

    T_reblock("an,bn,in,jn") = T_("a,b,i,j") * j_trans_array("j,jn")
                                     * i_trans_array("i,in")
                                     * b_trans_array("b,bn")
                                     * a_trans_array("a,an");



    std::cout << "first reblocking step worked" << std::endl;

    // Need to get number of tiles and number of elements
    // along each dim of T_reblock



    /// Steps (2) and (3): For each ti, tj pair, compute D_titj
    /// and transform to matrix
    /// Store D_titj matrices in vector





//    /// Steps (2) and (3): For each ti, tj pair, compute D_titj array and convert
//    // D_titj to local tensor

//    // Store local D_titj tensors in a vector with
//    // each D_titj tensor indexed by {ti,tj}
//    std::vector<Tile> D_titj_tensor_list(ntiles_i * ntiles_j);

//    typedef std::vector<std::size_t> block;

//    // Loop over ti, tj; get block indices
//    for (std::size_t ti = 0; ti != ntiles_i; ++ti) {
//      const std::size_t i_low = ti;
//      const std::size_t i_up = ti + 1;

//      for (std::size_t tj = 0; tj != ntiles_j; ++tj) {
//        const std::size_t j_low = tj;
//        const std::size_t j_up = tj + 1;

//        // auto delta_titj = (ti == tj) ? 1.0 : 0;

//        auto titj = ti * ntiles_j + tj; // ti,tj index for D_titj_tensor
//        std::cout << "Successfully evaluated titj" << std::endl;
      

//        // D_titj array
//        T D_titj_array;

//        // D_titj local tensor
//        Tile D_titj_tensor;

//    std::cout << "Successfully declared D_titj_array, D_titj_tensor" << std::endl;




//        // Lower and upper bounds for block expressions
//        block low_bound{0, 0, i_low, j_low};
//        block up_bound{ntiles_a, ntiles_b, i_up, j_up};

//        auto T_caij_block = T_("c,a,i,j").block(low_bound, up_bound);
//        auto T_caji_block = T_("c,a,j,i").block(low_bound, up_bound);
//        auto T_cbij_block = T_("c,b,i,j").block(low_bound, up_bound);
//        auto T_acij_block = T_("a,c,i,j").block(low_bound, up_bound);
//        auto T_acji_block = T_("a,c,j,i").block(low_bound, up_bound);
//        auto T_bcij_block = T_("b,c,i,j").block(low_bound, up_bound);

//        D_titj_array("a,b") =
//                                  ((4 * T_caij_block - 2 * T_caji_block) * T_cbij_block) +
//                                  ((4 * T_acij_block - 2 * T_acji_block) * T_bcij_block);


//        // D_titj_array("a,b,i,j") = T_caij_block;


//        //Block expression to form D from T
//        //** note ** divide by 1 + delta_ij later when D_ij matrix has been formed
//        // D_titj_array("a,b,i,j") =
//        //   (4 * T_("c,a,i,j").block(low_bound, up_bound) -
//        //    2 * T_("c,a,j,i").block(low_bound, up_bound)) *
//        //     T_("c,b,i,j").block(low_bound, up_bound) +
//        //   (4 * T_("a,c,i,j").block(low_bound, up_bound) -
//        //    2 * T_("a,c,j,i").block(low_bound, up_bound)) *
//        //     T_("b,c,i,j").block(low_bound, up_bound);

//      std::cout << "For ti = " << ti << ", tj = " << tj
//      << " block expression successful." << std::endl;

 



//        // // Divide D_titj_array by 1 + delta(ti,tj)
//        // if (ti == tj) {
//        //   D_titj_array = D_titj_array / (2.0);
//        // }

//        // Pack contents of D_titj_array into local tensor

//        // Loop over ta, tb
//        for (std::size_t ta = 0; ta != ntiles_a; ++ta) {

//          for (std::size_t tb = 0; tb != ntiles_b; ++tb) {

//            // Pull particular tile corresponding to ta, tb, ti, tj
//            Tile D_ab = D_titj_array.find({ta, tb, ti, tj});

//            // Element info along all four dims of D_ab
//            auto a0 = D_ab.range().lobound()[0];
//            auto an = D_ab.range().upbound()[0];
//            auto a_ext = an - a0;

//            auto b0 = D_ab.range().lobound()[1];
//            auto bn = D_ab.range().upbound()[1];
//            auto b_ext = bn - b0;

//            auto i0 = D_ab.range().lobound()[2];
//            auto in = D_ab.range().upbound()[2];
//            auto i_ext = in - i0;

//            auto j0 = D_ab.range().lobound()[3];
//            auto jn = D_ab.range().upbound()[3];
//            auto j_ext = jn - j0;

//            auto da = b_ext * i_ext * j_ext;
//            auto db = i_ext * j_ext;
//            auto di = j_ext;

//            // Element off-set in going from D_ab to D_titj tensor
//            auto a_off = ta * a_ext;
//            auto b_off = tb * b_ext;
//            auto i_off = ti * i_ext;
//            auto j_off = tj * j_ext;


//            // Loop over elements of D_ab
//            // a_shift is shifted a, etc.
//            for (auto a = a0; a != an; ++a) {
//              auto a_shift = a + a_off;

//              for (auto b = b0; b != bn; ++b) {
//                auto b_shift = b + b_off;

//                for (auto i = i0; i != in; ++i) {
//                  auto i_shift = i + i_off;

//                  for (auto j = j0; j != jn; ++j) {
//                    auto j_shift = j + j_off;

//                    auto D_ab_idx = a * da + b * db + i * di + j;
//                    auto D_titj_idx = a_shift * da + b_shift * db + i_shift * di + j_shift;

//                    D_titj_tensor[D_titj_idx] = D_ab[D_ab_idx];

//                  } // j
//                } // i
//              } // b
//            } // a
//          } // tb
//        } // ta

//        // Store D_titj_tensor in D_titj_tensor_list
//        D_titj_tensor_list[titj] = D_titj_tensor;

//      } // tj
//    } // ti

//    std::cout << "Successfully formed D_ij" << std::endl;



//    /// Step (4): For each {ti, tj} pair, for each {i,j} pair
//    // convert slice of D_ij to matrix, diagonalize, truncate PNOs,
//    // store matrix of PNOs, etc.

//    // Vector hold D_ij matrices for particular {i,j}
//    std::vector<Eigen::MatrixXd> D_ij_mat_list(nocc_act * nocc_act);

//    for (std::size_t ti = 0; ti != ntiles_i; ++ti) {
//      for (std::size_t tj = 0; tj != ntiles_j; ++tj) {
        
//        auto titj = ti * ntiles_j + tj;
//        Tile D_ij = D_titj_tensor_list[titj];

//        // Get number of a, b, i, j in D_ij

//        auto a0 = D_ij.range().lobound()[0];
//        auto an = D_ij.range().upbound()[0];
//        const std::size_t a_ext = an - a0;

//        auto b0 = D_ij.range().lobound()[1];
//        auto bn = D_ij.range().upbound()[1];
//        const std::size_t b_ext = bn - b0;

//        auto i0 = D_ij.range().lobound()[2];
//        auto in = D_ij.range().upbound()[2];
//        auto i_ext = in - i0;

//        auto j0 = D_ij.range().lobound()[3];
//        auto jn = D_ij.range().upbound()[3];
//        auto j_ext = jn - j0;

//        // Offset for index of D_ij
//        auto i_off = ti * i_ext;
//        auto j_off = tj * j_ext;

//        // Block indices for converting to matrix
//        auto a_low = a0;
//        auto a_hi = an;

//        auto b_low = b0;
//        auto b_hi = bn;

//        // Loop over individual i,j elements in D_titj tensor

//        for (auto i = i0; i != in; ++i) {
//          auto i_low = i;
//          auto i_hi = i + 1;
//          auto i_shift = i + i_off;

//          for (auto j = j0; j != jn; ++j) {
//            auto j_low = j;
//            auto j_hi = j + 1;
//            auto j_shift = j + j_off;

//            auto ij = i_shift * nocc_act + j_shift;
//            int delta_ij = (i == j) ? 1 : 0;

//            block low{a_low, b_low, i_low, j_low};
//            block high{a_hi, b_hi, i_hi, j_hi};

//            Tile D_ij_data = D_ij.block(low, high);

//            Eigen::MatrixXd D_ij_mat = TA::eigen_map(D_ij_data, a_ext, b_ext);
//            D_ij_mat = D_ij_mat / (1 + delta_ij);

//            D_ij_mat_list[ij] = D_ij_mat;

//          }
//        }
//      }
//    }


//    // Loop over i,j indices and diagonalize D_ij
//    for (auto i = 0; i != nocc_act; ++i) {
//      for (auto j = 0; j != nocc_act; ++j) {
//        auto ij = i * nocc_act + j;
//        Eigen::MatrixXd D_ij = D_ij_mat_list[ij];

//        // Diagonalize D_ij
//        es.compute(D_ij);
//        Eigen::MatrixXd pno_ij = es.eigenvectors();
//        auto occ_ij = es.eigenvalues();

//        // truncate PNOs
//        size_t pnodrop = 0;
//        if (tpno_ != 0.0) {
//          for (size_t k = 0; k != occ_ij.rows(); ++k) {
//            if (!(occ_ij(k) >= tpno_))
//              ++pnodrop;
//            else
//              break;
//          }
//        }
//        const auto npno = nuocc - pnodrop;

//        // Store truncated PNOs
//        Eigen::MatrixXd pno_trunc = pno_ij.block(0, pnodrop, nuocc, npno);
//        pnos_[ij] = pno_trunc;

//        // Transform F to PNO space
//        Eigen::MatrixXd F_pno_ij = pno_trunc.transpose() * F_uocc * pno_trunc;

//        // Store just the diagonal elements of F_pno_ij
//        F_pno_diag_[ij] = F_pno_ij.diagonal();

//        /////// Transform PNOs to canonical PNOs if pno_canonical_ == true

//        if (pno_canonical_ == "true" && npno > 0) {
//          // Compute eigenvectors of F in PNO space
//          es.compute(F_pno_ij);
//          Eigen::MatrixXd pno_transform_ij = es.eigenvectors();

//          // Transform pno_ij to canonical PNO space; pno_ij -> can_pno_ij
//          Eigen::MatrixXd can_pno_ij = pno_trunc * pno_transform_ij;

//          // Replace standard with canonical PNOs
//          pnos_[ij] = can_pno_ij;
//          F_pno_diag_[ij] = es.eigenvalues();
//        }

//        // truncate OSVs

//        // auto osvdrop = 0;
//        if (i == j) {
//          size_t osvdrop = 0;
//          if (tosv_ != 0.0) {
//            for (size_t k = 0; k != occ_ij.rows(); ++k) {
//              if (!(occ_ij(k) >= tosv_))
//                ++osvdrop;
//              else
//                break;
//            }
//          }
//          const auto nosv = nuocc - osvdrop;
//          if (nosv == 0) {  // all OSV truncated indicates total nonsense
//            throw LimitExceeded<size_t>("all OSVs truncated", __FILE__,
//                                        __LINE__, 1, 0);
//          }

//          // Store truncated OSVs
//          Eigen::MatrixXd osv_trunc = pno_ij.block(0, osvdrop, nuocc, nosv);
//          osvs_[i] = osv_trunc;

//          // Transform F to OSV space
//          Eigen::MatrixXd F_osv_i = osv_trunc.transpose() * F_uocc * osv_trunc;

//          // Store just the diagonal elements of F_osv_i
//          F_osv_diag_[i] = F_osv_i.diagonal();

//          /////// Transform OSVs to canonical OSVs if pno_canonical_ == true
//          if (pno_canonical_ == "true") {
//            // Compute eigenvectors of F in OSV space
//            es.compute(F_osv_i);
//            Eigen::MatrixXd osv_transform_i = es.eigenvectors();

//            // Transform osv_i to canonical OSV space: osv_i -> can_osv_i
//            Eigen::MatrixXd can_osv_i = osv_trunc * osv_transform_i;

//            // Replace standard with canonical OSVs
//            osvs_[i] = can_osv_i;
//            F_osv_diag_[i] = es.eigenvalues();
//          }
//        } // if (i == j)
//      } // j
//    } // i

//    auto sum_osv = 0;
//    for (int i = 0; i < nocc_act; ++i) {
//      sum_osv += osvs_[i].cols();
//    }
//    auto ave_nosv = sum_osv / nocc_act;
//    ExEnv::out0() << "The average number of OSVs is " << ave_nosv << std::endl;

//    auto sum_pno = 0;
//    for (int i = 0; i < nocc_act; ++i) {
//      for (int j = 0; j < nocc_act; ++j) {
//        sum_pno += pnos_[i * nocc_act + j].cols();
//      }
//    }
//    auto ave_npno = sum_pno / (nocc_act * nocc_act);
//    ExEnv::out0() << "The average number of PNOs is " << ave_npno << std::endl;
//  }





  



 /// !!! Original PNO formation code !!! ///
 // Do not delete! //


     // zero out amplitudes
     if (!T_.is_initialized()) {
       T_ = Array(world, K.trange(), K.shape());
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
     const auto C_i_eig = TA::array_to_eigen(ofac.retrieve("i").coefs());
     const auto C_a_eig = TA::array_to_eigen(ofac.retrieve("a").coefs());
     const auto libint2_shells = fac.basis_registry()->retrieve(L"μ")->flattened_shells();

     // write out active occupied orbitals
     auto occs = Eigen::VectorXd::Constant(C_i_eig.cols(), 2.0);
     libint2::molden::Export xport(libint2_atoms, libint2_shells, C_i_eig, occs);
     xport.write("occ.molden");
 #endif

  

    


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

         // Store truncated PNOs
         // pnos[i*nocc_act + j] = pno_ij.block(0,pnodrop,nuocc,npno);
         Eigen::MatrixXd pno_trunc = pno_ij.block(0, pnodrop, nuocc, npno);
         // pnos_[ij] = pno_trunc;
         pnos_[i * nocc_act + j] = pno_trunc;

 #if PRODUCE_PNO_MOLDEN_FILES
         // write PNOs to Molden
         {
           Eigen::MatrixXd molden_coefs(C_i_eig.rows(), 2 + pno_trunc.cols());
           molden_coefs.col(0) = C_i_eig.col(i);
           molden_coefs.col(1) = C_i_eig.col(j);
           molden_coefs.block(0, 2, C_i_eig.rows(), pno_trunc.cols()) = C_a_eig * pno_trunc;

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

         if (pno_canonical_ == "true" && npno > 0) {
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

         // auto nosv = 0;

         // auto osvdrop = 0;
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
           Eigen::MatrixXd F_osv_i = osv_trunc.transpose() * F_uocc * osv_trunc;

           // Store just the diagonal elements of F_osv_i
           F_osv_diag_[i] = F_osv_i.diagonal();

           /////// Transform OSVs to canonical OSVs if pno_canonical_ == true
           if (pno_canonical_ == "true") {
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
     ExEnv::out0() << "The average number of OSVs is " << ave_nosv << std::endl;

     auto sum_pno = 0;
     for (int i = 0; i < nocc_act; ++i) {
       for (int j = 0; j < nocc_act; ++j) {
         sum_pno += pnos_[i * nocc_act + j].cols();
       }
     }
     auto ave_npno = sum_pno / (nocc_act * nocc_act);
     ExEnv::out0() << "The average number of PNOs is " << ave_npno << std::endl;
   }



  virtual ~PNOSolver() = default;

  /// @return PNO truncation threshold
  double tpno() const { return tpno_; }
  /// @return OSV truncation threshold
  double tosv() const { return tosv_; }

  const auto& pno(int i, int j) const { return pnos_[i * nocc_act_ + j]; }
  const auto& osv(int i) const { return osvs_[i]; }

 private:
  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space
  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    auto delta_t1_ai = jacobi_update_t1(r1, F_occ_act_, F_osv_diag_, osvs_);
    auto delta_t2_abij = jacobi_update_t2(r2, F_occ_act_, F_pno_diag_, pnos_);
    t1("a,i") += delta_t1_ai("a,i");
    t2("a,b,i,j") += delta_t2_abij("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }

  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    update_only(t1, t2, r1, r2);
    T r1_osv = osv_transform_ai(r1, osvs_);
    T r2_pno = pno_transform_abij(r2, pnos_);
    mpqc::cc::T1T2<T, T> r(r1_osv, r2_pno);
    mpqc::cc::T1T2<T, T> t(t1, t2);
    this->diis().extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
  }





  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t2(
      const TA::DistArray<Tile, Policy>& r2_abij,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_pno_diag,
      const std::vector<Eigen::MatrixXd>& pnos) {

    auto update2 = [F_occ_act, F_pno_diag, pnos, this](
                       Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      auto ij = i * nocc_act_ + j;
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
      Eigen::MatrixXd delta_t2_full =
          pno_ij * delta_t2_pno * pno_ij.transpose();

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

    auto delta_t2_abij = TA::foreach(r2_abij, update2);
    delta_t2_abij.world().gop.fence();
    return delta_t2_abij;
  }



  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t1(
      const TA::DistArray<Tile, Policy>& r1_ai,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_osv_diag,
      const std::vector<Eigen::MatrixXd>& osvs) {
    auto update1 = [F_occ_act, F_osv_diag, osvs](Tile& result_tile,
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
      // basis this matrix will then be back transformed to full basis before
      // being converted to a tile
      Eigen::VectorXd delta_t1_osv = r1_osv;

      // Select correct vector containing diagonal elements of Fock matrix in
      // OSV basis
      const Eigen::VectorXd& ens_uocc = F_osv_diag[i];

      // Determine number of OSVs
      const auto nosv = ens_uocc.rows();

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
      // Eigen::MatrixXd delta_t1_full = osv_i * delta_t1_osv *
      // osv_i.transpose();
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

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> pno_transform_abij(
      const TA::DistArray<Tile, Policy>& abij,
      const std::vector<Eigen::MatrixXd>& pnos) {

    auto tform = [pnos, this](
        Tile& result_tile, const Tile& arg_tile) {

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      const auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd pno_ij = pnos[ij];
      const auto nuocc = pno_ij.rows();
      const auto npno = pno_ij.cols();

      // Convert data in tile to Eigen::Map and transform to PNO basis
      const Eigen::MatrixXd result_eig =
          pno_ij.transpose() * TA::eigen_map(arg_tile, nuocc, nuocc) * pno_ij;

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{npno,npno,1l,1l});
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

    auto result = TA::foreach(abij, tform);
    result.world().gop.fence();
    return result;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> osv_transform_ai(
      const TA::DistArray<Tile, Policy>& ai,
      const std::vector<Eigen::MatrixXd>& osvs) {

    auto tform = [osvs, this](
        Tile& result_tile, const Tile& arg_tile) {

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
      result_tile = Tile(TA::Range{nosv,1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nosv; ++r) {
        const auto elem = result_eig(r, 0);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach(ai, tform);
    result.world().gop.fence();
    return result;
  }

  // squared norm of 1-body residual in OSV subspace
  struct R1SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R1SquaredNormReductionOp(PNOSolver<T>* solver) : solver_(solver) {}

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

    PNOSolver<T>* solver_;
  };  // R1SquaredNormReductionOp

  // squared norm of 2-body residual in PNO subspace
  struct R2SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R2SquaredNormReductionOp(PNOSolver<T>* solver) : solver_(solver) {}

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

    PNOSolver<T>* solver_;
  };  // R2SquaredNormReductionOp

 public:
  /// Overrides Solver<T,T>::error()
  virtual double error(const T& r1, const T& r2) override {
    R1SquaredNormReductionOp op1(this);
    R2SquaredNormReductionOp op2(this);
    return sqrt(r1("a,i").reduce(op1).get() + r2("a,b,i,j").reduce(op2).get()) /
           (size(r1) + size(r2));
  }

 private:
  Factory<T>& factory_;
  std::string pno_method_;     //!< the PNO construction method
  std::string pno_canonical_;  //!< whether or not to canonicalize PNO/OSV
  double tpno_;                //!< the PNO truncation threshold
  double tosv_;                //!< the OSV (diagonal PNO) truncation threshold
  int nocc_act_;               //!< the number of active occupied orbitals
  Array T_;                    //!< the array of MP2 T amplitudes
  Array D_;                    //!< the array of MP2 density values 

  Eigen::MatrixXd F_occ_act_;

  // For storing PNOs and and the Fock matrix in the PNO basis
  std::vector<Eigen::MatrixXd> pnos_;
  std::vector<Eigen::VectorXd> F_pno_diag_;

  // For storing OSVs (PNOs when i = j) and the Fock matrix in
  // the OSV basis
  std::vector<Eigen::MatrixXd> osvs_;
  std::vector<Eigen::VectorXd> F_osv_diag_;
};  // class: PNO solver


////// SVO2 solver ////////

/// SVOSolver updates the CC T amplitudes using standard Jacobi+DIIS in SVO2
/// space
/// @warning This class assumes that the 1- and 2-body amplitudes/residuals
///          given to Solver::update() are laid out as "a,i" and "a,b,i,j",
///          respectively
template <typename T>
class SVOSolver : public ::mpqc::cc::DIISSolver<T, T>,
                  public madness::WorldObject<SVOSolver<T>> {
 public:
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
  SVOSolver(const KeyVal& kv, Factory<T>& factory)
      : ::mpqc::cc::DIISSolver<T, T>(kv),
        madness::WorldObject<SVOSolver<T>>(factory.world()),
        factory_(factory),
        tsvo2_(kv.value<double>("tsvo2", 1.e-5)),
        tsvo1_(kv.value<double>("tsvo1", 1.e-7)) {
    // part of WorldObject initialization
    this->process_pending();

    // Check that tiling is done appropriately
    if (kv.exists("occ_block_size")) {
      int occ_block_size_ = (kv.value<int>("occ_block_size"));
      if (occ_block_size_ != 1) {
        throw InputError("occ_block_size must be set to 1 in the input file.");
      }
    } 
    else {
      throw InputError("occ_block_size was not specified in the input file.");
    }

    if (kv.exists("unocc_block_size")) {
      int unocc_block_size_ = (kv.value<int>("unocc_block_size"));
      if (unocc_block_size_ < 1000000000) {
        throw InputError(
            "unocc_block_size must be greater than or equal to 1000000000 in "
            "the input file.");
      }
    } 
    else {
      throw InputError("unocc_block_size was not specified in the input file.");
    }

    // Compute integrals

    auto& fac = factory_;
    auto& world = fac.world();
    auto& ofac = fac.orbital_registry();

    auto nocc = ofac.retrieve("m").rank();
    auto nocc_act = ofac.retrieve("i").rank();
    nocc_act_ = nocc_act;
    auto nuocc = ofac.retrieve("a").rank();
    auto nfzc = nocc - nocc_act;

    // Form Fock array
    auto F = fac.compute(L"<p|F|q>[df]");

    // Select just diagonal elements of Fock aray and transform
    // to Eigen vector; use for computing SVO2s
    Eigen::VectorXd eps_p = TA::array_to_eigen(F).diagonal();
    auto eps_o = eps_p.segment(nfzc, nocc_act);
    auto eps_v = eps_p.tail(nuocc);

    // Transform entire Fock array to Eigen Matrix
    Eigen::MatrixXd F_all = TA::array_to_eigen(F);

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
      T_ = Array(world, K.trange(), K.shape());
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

          } // for each b
        } // for each a

        // Compute SVD of each T_ij matrix
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T_ij, Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto sing_vals = svd.singularValues();

        // truncate SVO2s
        size_t svo2_keep = 0;
        if (tsvo2_ != 0.0) {
          for (size_t k = 0; k != sing_vals.rows(); ++k) {
            if (sing_vals(k) >= tsvo2_)
              ++svo2_keep;
            else
              break;
          } // for each k
          if (svo2_keep == 0) {  // all SVO2 truncated indicates total nonsense
            throw LimitExceeded<size_t>("all SVO2s truncated", __FILE__,
                                        __LINE__, 1, 0);
          }
        } // if tsvo2 != 0

        else {  // handles case where tsvo2_ == 0
          svo2_keep = nuocc;
        }

        const auto nsvo2 = svo2_keep;
        //std::cout << "nsvo2 = " << nsvo2 << std::endl;

        // store truncated SVO2s
        Eigen::MatrixXd r_svo2_trunc = svd.matrixV().block(0, 0, nuocc, nsvo2);
        Eigen::MatrixXd l_svo2_trunc = svd.matrixU().block(0, 0, nuocc, nsvo2);
        r_svo2s_[i*nocc_act + j] = r_svo2_trunc;
        l_svo2s_[i*nocc_act + j] = l_svo2_trunc;

        // transform F to right SVO2 space and store just diagonal elements
        //Eigen::MatrixXd F_r_svo2 = svd.matrixV().transpose() * F_uocc * svd.matrixV();
        Eigen::MatrixXd F_r_svo2 = r_svo2_trunc.transpose() * F_uocc * r_svo2_trunc;
        F_r_svo2_diag_[i*nocc_act + j] = F_r_svo2.diagonal();

        // transform F to left SVO2 space and store just diagonal elements
        //Eigen::MatrixXd F_l_svo2 = svd.matrixU().transpose() * F_uocc.transpose() * svd.matrixU();
        Eigen::MatrixXd F_l_svo2 = l_svo2_trunc.transpose() * F_uocc.transpose() * l_svo2_trunc;
        F_l_svo2_diag_[i*nocc_act + j] = F_l_svo2.diagonal();



        // if i == j, truncate and store SVO1s
        if (i == j) {
          size_t svo1_keep = 0;
          if (tsvo1_ != 0.0) {
            for (size_t k = 0; k != sing_vals.rows(); ++k) {
              if (sing_vals(k) >= tsvo1_)
                ++svo1_keep;
              else
                break;
            } // for each k
            if (svo1_keep == 0) {  // all SVO1 truncated indicates total nonsense
              throw LimitExceeded<size_t>("all SVO1s truncated", __FILE__,
                                          __LINE__, 1, 0);
            }
          } // if tsvo1 != 0

          else {  // handles case where tsvo1_ == 0
            svo1_keep = nuocc;
          }

          const auto nsvo1 = svo1_keep;
          //std::cout << "nsvo1 = " << nsvo1 << std::endl;

          // store truncated SVO1s
          Eigen::MatrixXd svo1_trunc = svd.matrixV().block(0, 0, nuocc, nsvo1);
          svo1s_[i] = svo1_trunc;

          // transform F to SVO1 space and store just diagonal elements
          Eigen::MatrixXd F_svo1 = svo1_trunc.transpose() * F_uocc * svo1_trunc;
          F_svo1_diag_[i] = F_svo1.diagonal();

        } // if (i == j)


      } // for each j
    } // for each i

    // Compute average number of SVO2s per pair and print out

    auto sum_svo2 = 0;
    for (int i=0; i<nocc_act; ++i) {
      for (int j=0; j<nocc_act; ++j) {
        sum_svo2 += r_svo2s_[i*nocc_act + j].cols();
      }
    } 
    auto ave_nsvo2 = sum_svo2 / (nocc_act * nocc_act);
    ExEnv::out0() << "The average number of SVO2s is " << ave_nsvo2 << std::endl;


    // Compute average number of SVO1s per pair and print out

    auto sum_svo1 = 0;
    for (int i=0; i<nocc_act; ++i) {
      sum_svo1 += svo1s_[i].cols();
    } 
    auto ave_nsvo1 = sum_svo1 / nocc_act;
    //std::cout << "sum_svo1 = " << sum_svo1 << std::endl;
    ExEnv::out0() << "The average number of SVO1s is " << ave_nsvo1 << std::endl;


  } // SVOSolver

  virtual ~SVOSolver() = default;

  /// @return SVO2 truncation threshold
  double tsvo2() const { return tsvo2_; }
  /// @return SVO1 truncation threshold
  double tsvo1() const { return tsvo1_; }

  const auto& l_svo2(int i, int j) const { return l_svo2s_[i*nocc_act_ + j]; }
  const auto& r_svo2(int i, int j) const { return r_svo2s_[i*nocc_act_ + j]; }

private:
  /// Overrides DIISSolver::update_only() .
  /// @note must override DIISSolver::update() also since the update must be
  ///      followed by backtransform updated amplitudes to the full space
  void update_only(T& t1, T& t2, const T& r1, const T& r2) override {
    auto delta_t1_ai = jacobi_update_t1(r1, F_occ_act_, F_svo1_diag_, svo1s_);

    auto delta_t2_abij = jacobi_update_t2(r2, F_occ_act_, F_l_svo2_diag_,
                                          F_r_svo2_diag_, l_svo2s_, r_svo2s_);
    t1("a,i") += delta_t1_ai("a,i");
    t2("a,b,i,j") += delta_t2_abij("a,b,i,j");
    t1.truncate();
    t2.truncate();
  }

  void update(T& t1, T& t2, const T& r1, const T& r2) override {
    update_only(t1, t2, r1, r2);
    T r1_svo1 = svo_transform_ai(r1, svo1s_);
    T r2_svo2 = svo_transform_abij(r2, l_svo2s_, r_svo2s_);
    mpqc::cc::T1T2<T, T> r(r1_svo1, r2_svo2);
    mpqc::cc::T1T2<T, T> t(t1, t2);
    this->diis().extrapolate(t, r);
    t1 = t.t1;
    t2 = t.t2;
  }

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t2(
      const TA::DistArray<Tile, Policy>& r2_abij,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_l_svo2_diag,
      const std::vector<Eigen::VectorXd>& F_r_svo2_diag,
      const std::vector<Eigen::MatrixXd>& l_svo2s,
      const std::vector<Eigen::MatrixXd>& r_svo2s) {

    auto update2 = [F_occ_act, F_l_svo2_diag, F_r_svo2_diag, l_svo2s, r_svo2s, this](
                       Tile& result_tile, const Tile& arg_tile) {

      result_tile = Tile(arg_tile.range());

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of PNOs
      auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd l_svo2_ij = l_svo2s[ij];
      Eigen::MatrixXd svo1_ij = r_svo2s[ij];

      // Extent data of tile
      const auto ext = arg_tile.range().extent_data();

      // Convert data in tile to Eigen::Map and transform to SVO2 basis
      const Eigen::MatrixXd r2_svo2 =
          l_svo2_ij.transpose() *
          TA::eigen_map(arg_tile, ext[0] * ext[2], ext[1] * ext[3]) * svo1_ij;

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
          l_svo2_ij * delta_t2_svo2 * svo1_ij.transpose();

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

    auto delta_t2_abij = TA::foreach(r2_abij, update2);
    delta_t2_abij.world().gop.fence();
    return delta_t2_abij;
  } // jacobi_update_t2

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> jacobi_update_t1(
      const TA::DistArray<Tile, Policy>& r1_ai,
      const Eigen::MatrixXd& F_occ_act,
      const std::vector<Eigen::VectorXd>& F_svo1_diag,
      const std::vector<Eigen::MatrixXd>& svo1s) {
    auto update1 = [F_occ_act, F_svo1_diag, svo1s, this](
                      Tile& result_tile, const Tile& arg_tile) {

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
      //Eigen::VectorXd delta_t1_full = delta_t1_svo1 * svo1_i;

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
  } // jacobi_update_t1

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> svo_transform_abij(
      const TA::DistArray<Tile, Policy>& abij,
      const std::vector<Eigen::MatrixXd>& l_svo2s,
      const std::vector<Eigen::MatrixXd>& r_svo2s) {

    auto tform = [l_svo2s, r_svo2s, this](
        Tile& result_tile, const Tile& arg_tile) {

      // determine i and j indices
      const auto i = arg_tile.range().lobound()[2];
      const auto j = arg_tile.range().lobound()[3];

      // Select appropriate matrix of SVO2s
      const auto ij = i * nocc_act_ + j;
      Eigen::MatrixXd l_svo2_ij = l_svo2s[ij];
      Eigen::MatrixXd svo1_ij = r_svo2s[ij];
      const auto nuocc = l_svo2_ij.rows();
      const auto nsvo2 = l_svo2_ij.cols();

      // Convert data in tile to Eigen::Map and transform to SVO2 basis
      const Eigen::MatrixXd result_eig =
          l_svo2_ij.transpose() * TA::eigen_map(arg_tile, nuocc, nuocc) * svo1_ij;

      // Convert result_eig to tile and compute norm
      result_tile = Tile(TA::Range{nsvo2,nsvo2,1l,1l});
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

    auto result = TA::foreach(abij, tform);
    result.world().gop.fence();
    return result;
  } // svo_transform_abij

  template <typename Tile, typename Policy>
  TA::DistArray<Tile, Policy> svo_transform_ai(
      const TA::DistArray<Tile, Policy>& ai,
      const std::vector<Eigen::MatrixXd>& svo1s) {

    auto tform = [svo1s, this](
        Tile& result_tile, const Tile& arg_tile) {

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
      result_tile = Tile(TA::Range{nsvo1,1l});
      typename Tile::scalar_type norm = 0.0;
      for (auto r = 0; r < nsvo1; ++r) {
        const auto elem = result_eig(r, 0);
        const auto abs_elem = std::abs(elem);
        norm += abs_elem * abs_elem;
        result_tile[r] = elem;
      }

      return std::sqrt(norm);
    };

    auto result = TA::foreach(ai, tform);
    result.world().gop.fence();
    return result;
  } // svo_transform_ai

  // squared norm of 1-body residual in SVO1 subspace
  struct R1SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R1SquaredNormReductionOp(SVOSolver<T>* solver) : solver_(solver) {}

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

    SVOSolver<T>* solver_;
  };  // R1SquaredNormReductionOp

  // squared norm of 2-body residual in SVO2 subspace
  struct R2SquaredNormReductionOp {
    // typedefs
    typedef typename TA::detail::scalar_type<T>::type result_type;
    typedef typename T::value_type argument_type;

    R2SquaredNormReductionOp(SVOSolver<T>* solver) : solver_(solver) {}

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

    SVOSolver<T>* solver_;
  };  // R2SquaredNormReductionOp

 public:
  /// Overrides Solver<T,T>::error()
  virtual double error(const T& r1, const T& r2) override {
    R1SquaredNormReductionOp op1(this);
    R2SquaredNormReductionOp op2(this);
    return sqrt(r1("a,i").reduce(op1).get() + r2("a,b,i,j").reduce(op2).get()) /
           (size(r1) + size(r2));
  }

  private:
  Factory<T>& factory_;
  double tsvo2_;        //!< the truncation threshold for SVO2s
  double tsvo1_;        //!< the truncation threshold for SVO1s
  int nocc_act_;        //!< the number of active occupied orbitals
  Array T_;

  Eigen::MatrixXd F_occ_act_;

  // For storing SVO2s and the Fock matrix in the SVO2 basis
  std::vector<Eigen::MatrixXd> l_svo2s_;
  std::vector<Eigen::MatrixXd> r_svo2s_;
  std::vector<Eigen::VectorXd> F_l_svo2_diag_;
  std::vector<Eigen::VectorXd> F_r_svo2_diag_;

  // For storing SVO1s and the Fock matrix in the SVO1 basis
  std::vector<Eigen::MatrixXd> svo1s_;
  std::vector<Eigen::VectorXd> F_svo1_diag_;



}; // class: SVOSolver



}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_LCAO_CC_SOLVERS_H_ */
