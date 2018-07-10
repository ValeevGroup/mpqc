//
// Created by Chong Peng on 8/21/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_T_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_T_H_

#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/misc/print.h"

// Eigen library
#include "mpqc/math/external/eigen/eigen.h"

// Laplace transformation of 2-electron integrals and amplitudes (t2 & t1)
#include "mpqc/chemistry/qc/lcao/cc/laplace_transform.h"

#include "mpqc/math/quadrature/gaussian.h"

namespace mpqc {
namespace lcao {

//

/**
 *  \brief CCSD_T class that compute CCSD(T) triple calculation
 *
 *  keyword to call this class CCSD(T)
 *
 */

template <typename Tile, typename Policy>
class CCSD_T : virtual public CCSD<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;

 private:
  /// if reblock occ and unocc space
  bool reblock_;

  /// if reblock inner contraction occ and unocc space
  bool reblock_inner_;

  /// if replicate integral g_cijk
  bool replicate_ijka_;

  /// string represent (T) approach, see keyval constructor for detail
  std::string approach_;

  /// occ reblock size
  std::size_t occ_block_size_;

  /// unocc reblock size
  std::size_t unocc_block_size_;

  /// inner contraction space reblock size
  std::size_t inner_block_size_;

  /// inner occ space TiledRange1
  TA::TiledRange1 tr_occ_inner_;

  /// inner unocc space TiledRange1
  TA::TiledRange1 tr_vir_inner_;

  /// increase size in unocc space in iteration of fine grain approach
  std::size_t increase_;

  /// num
  std::size_t n_laplace_quad_;

  /// (T) energy
  double triples_energy_;

  /// local TRange1Engine to allow reblocking
  std::shared_ptr<const ::mpqc::utility::TRange1Engine> trange1_engine_;
  const std::shared_ptr<const ::mpqc::utility::TRange1Engine>
      &local_trange1_engine() const {
    return trange1_engine_;
  }

 public:
  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords : all keywords for CCSD
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | @c approach | string | @c coarse | the (T) algorithm; valid choices are <ul> <li>@c coarse (parallelized over {a,b,c}, each assigned round-robin to single node) <li/> @c fine (same as @coarse, but each {a,b,c} task is executed over the entire machine; less scalable than @coarse ) <li/> @c straight (for reference purposes only) <li/> @c laplace (Laplace MO-based implementation) </ul> |
   * | @c increase | int | @c 2 | number of unoccupied tiles (per dimension) to load in each abc loop; valid only if @c approach=fine  |
   * | @c reblock_occ | int | @c 8 | the block size used for the occupied orbitals |
   * | @c reblock_unocc | int | @c 8 | the block size used for the unoccupied orbitals |
   * | @c reblock_inner | int | number of orbitals | the block size for the inner (contraction) dimension; set to 0 to disable reblock inner; only used if @c approach=laplace |
   * | @c replicate_ijka | bool | @c false | whether to replicate integral <ij\|ka>, the smallest 2-body integral in (T); valid only with @c approach=coarse |
   * | @c quadrature_points | int | @c 4 | number of quadrature points for the Laplace transform; valid only if @c approach=laplace |
   */
  // clang-format on

  CCSD_T(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    reblock_ = kv.exists("reblock_occ") || kv.exists("reblock_unocc");

    occ_block_size_ = kv.value<int>("reblock_occ", 8);
    unocc_block_size_ = kv.value<int>("reblock_unocc", 8);
    replicate_ijka_ = kv.value<bool>("replicate_ijka", false);

    // default value is size of total number of orbitals, which makes it 1 block
    std::size_t n_obs = this->wfn_world()
                            ->basis_registry()
                            ->retrieve(OrbitalIndex(L"κ"))
                            ->nfunctions();
    inner_block_size_ = reblock_ ? kv.value<int>("reblock_inner", n_obs) : 0;
    reblock_inner_ = (inner_block_size_ == 0) ? false : true;

    increase_ = kv.value<int>("increase", 2);
    approach_ = kv.value<std::string>("approach", "coarse");
    if (approach_ != "coarse" && approach_ != "fine" &&
        approach_ != "straight" && approach_ != "laplace") {
      throw InputError("Invalid (T) approach! \n", __FILE__, __LINE__,
                       "approach");
    }

    // no default reblock inner if use laplace
    if (approach_ == "laplace") {
      reblock_ = false;
      reblock_inner_ = false;
    }

    // set number of quadrature points for the Laplace transform
    n_laplace_quad_ = kv.value<int>("quadrature_points", 4);
  }

  virtual ~CCSD_T() {}

  void obsolete() override {
    triples_energy_ = 0.0;
    CCSD<Tile, Policy>::obsolete();
  }

  double triples_energy() const { return this->triples_energy_; }

 protected:
  void compute_ccsd_t() {
    auto &world = this->wfn_world()->world();
    auto time0 = mpqc::fenced_now(world);

    // clean all LCAO integral
    this->lcao_factory().registry().purge();

    if (approach_ != "laplace") {
      // reblock occ and unocc space
      if (reblock_ || reblock_inner_) {
        reblock();
      }
    }

    // start CCSD(T)
    ExEnv::out0() << "\nBegining CCSD(T) " << std::endl;
    TArray t1 = this->t1();
    TArray t2 = this->t2();
    if (approach_ == "coarse") {
      triples_energy_ = compute_ccsd_t_coarse_grain(t1, t2);
    } else if (approach_ == "fine") {
      triples_energy_ = compute_ccsd_t_fine_grain(t1, t2);
    } else if (approach_ == "straight") {
      triples_energy_ = compute_ccsd_t_straight(t1, t2);
    } else if (approach_ == "laplace") {
      triples_energy_ = compute_ccsd_t_laplace_transform(t1, t2);
    }

    auto time1 = mpqc::fenced_now(world);
    auto duration1 = mpqc::duration_in_s(time0, time1);

    ExEnv::out0() << "(T) Energy: " << triples_energy_ << " Time: " << duration1
                  << " S \n";
  }

  void evaluate(Energy *result) override {
    if (!this->computed()) {
      auto &world = this->lcao_factory().world();

      CCSD<Tile, Policy>::evaluate(result);
      double ccsd_energy = this->get_value(result).derivs(0)[0];

      auto time0 = mpqc::fenced_now(world);
      // compute
      compute_ccsd_t();

      this->computed_ = true;
      this->set_value(result, ccsd_energy + triples_energy_);
      auto time1 = mpqc::fenced_now(world);
      auto duration0 = mpqc::duration_in_s(time0, time1);
      ExEnv::out0() << "(T) Time in CCSD(T): " << duration0 << " S"
                    << std::endl;
    }
  }

 private:
  double compute_ccsd_t_coarse_grain(TArray &t1, TArray &t2) {
    auto &global_world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    // get integral
    TArray g_aijk = get_aijk();
    TArray g_abij = get_abij();
    TArray g_dabi = get_abci();

    // T2
    TArray t2_left = t2;
    TArray t2_right = t2;

    // if reblock_inner, needs two copy of T2 with different blocking
    if (reblock_inner_) {
      reblock_inner_t2(t2_left, t2_right);
    }

    auto trange1_engine = this->local_trange1_engine() == nullptr
                          ? this->trange1_engine()
                          : this->local_trange1_engine();
    // get trange1
    auto tr_occ = trange1_engine->get_active_occ_tr1();
    auto tr_vir = trange1_engine->get_vir_tr1();

    auto n_tr_occ = trange1_engine->get_active_occ_blocks();
    auto n_tr_vir = trange1_engine->get_vir_blocks();

    // TiledRange1 for occ, unocc and inner contraction space
    auto n_tr_occ_inner = n_tr_occ;
    auto n_tr_vir_inner = n_tr_vir;
    if (reblock_inner_) {
      n_tr_occ_inner = tr_occ_inner_.tiles_range().second;
      n_tr_vir_inner = tr_vir_inner_.tiles_range().second;
    }

    std::size_t vir_block_size =
        trange1_engine->get_vir_block_size();
    std::size_t n_occ = trange1_engine->get_active_occ();
            std::size_t n_blocks = n_tr_occ * n_tr_occ * n_tr_occ;
            double mem = (n_occ * n_occ * n_occ * vir_block_size * vir_block_size *
                vir_block_size * 8) /
                1.0e9;

            ExEnv::out0() << "Number of blocks at each iteration: " << n_blocks
                          << std::endl;
            ExEnv::out0() << "Size of T3 or V3 at each iteration per node: 3x" << mem
                          << " GB" << std::endl;

            // split global_world
            const auto rank = global_world.rank();
            const auto size = global_world.size();

            madness::World *tmp_ptr;
            std::shared_ptr<madness::World> world_ptr;

            if (size > 1) {
              SafeMPI::Group group = global_world.mpi.comm().Get_group().Incl(1, &rank);
              SafeMPI::Intracomm comm = global_world.mpi.comm().Create(group);
              world_ptr = std::make_shared<madness::World>(comm);
              tmp_ptr = world_ptr.get();
            } else {
              tmp_ptr = &global_world;
            }
            auto &this_world = *tmp_ptr;
            global_world.gop.fence();

            TArray g_cjkl_global(
                this_world, g_aijk.trange(), g_aijk.shape(),
                Policy::default_pmap(this_world,
                                     g_aijk.trange().tiles_range().volume()));

            TArray t1_this(
                this_world, t1.trange(), t1.shape(),
                Policy::default_pmap(this_world, t1.trange().tiles_range().volume()));

            // by default replicate t1 array
            if (size > 1) {
              // replicate t1 array
              t1_this("a,i") = t1("a,i").set_world(this_world);
            } else {
              // on 1 node, no need to replicate t1
              t1_this = t1;
            }

            // based on replicate keyword, replicate g_cjkl array
            if (replicate_ijka_ && size > 1) {
              // replicate g_cjkl
              g_cjkl_global("c,j,k,l") = g_aijk("c,j,k,l").set_world(this_world);
            } else {
              // don't replicate g_cjkl
              g_cjkl_global = g_aijk;
            }

            typedef std::vector<std::size_t> block;
            // lambda function to compute t3

            auto compute_t3_c = [&g_dabi, &t2_right, &t2_left, &g_cjkl_global,
                n_tr_vir_inner, n_tr_occ,
                n_tr_occ_inner](std::size_t a, std::size_t b,
                                std::size_t c, TArray &t3) {
              // index
              std::size_t a_low = a;
              std::size_t a_up = a + 1;
              std::size_t b_low = b;
              std::size_t b_up = b + 1;
              std::size_t c_low = c;
              std::size_t c_up = c + 1;

              {
                // block for g_dbai
                block g_dabi_low{0, a_low, b_low, 0};
                block g_dabi_up{n_tr_vir_inner, a_up, b_up, n_tr_occ};

                auto block_g_dabi = g_dabi("d,a,b,i").block(g_dabi_low, g_dabi_up);

                // block for t2_dcjk
                block t2_dcjk_low{0, c_low, 0, 0};
                block t2_dcjk_up{n_tr_vir_inner, c_up, n_tr_occ, n_tr_occ};
                auto block_t2_dcjk = t2_left("d,c,j,k").block(t2_dcjk_low, t2_dcjk_up);

                // block for g_cjkl
                block g_cjkl_low{c_low, 0, 0, 0};
                block g_cjkl_up{c_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};
                TArray g_cjkl;
                auto block_g_cjkl =
                    g_cjkl_global("c,j,k,l").block(g_cjkl_low, g_cjkl_up);

                // block for t2_abil
                block t2_abil_low{a_low, b_low, 0, 0};
                block t2_abil_up{a_up, b_up, n_tr_occ, n_tr_occ_inner};

                auto block_t2_abil = t2_right("a,b,i,l").block(t2_abil_low, t2_abil_up);

                if (t3.is_initialized()) {
                  t3("a,b,i,c,j,k") +=
                      (block_g_dabi * block_t2_dcjk) - (block_t2_abil * block_g_cjkl);
                } else {
                  t3("a,b,i,c,j,k") =
                      (block_g_dabi * block_t2_dcjk) - (block_t2_abil * block_g_cjkl);
                }
              }
            };

            auto compute_t3 = [&g_dabi, &t2_right, n_tr_vir_inner, n_tr_occ,
                n_tr_occ_inner](std::size_t a, std::size_t b, TArray &t3,
                                const TArray &this_t2_dcjk,
                                const TArray &this_g_cjkl) {
              // index
              std::size_t a_low = a;
              std::size_t a_up = a + 1;
              std::size_t b_low = b;
              std::size_t b_up = b + 1;

              {
                // block for g_dbai
                block g_dabi_low{0, a_low, b_low, 0};
                block g_dabi_up{n_tr_vir_inner, a_up, b_up, n_tr_occ};

                auto block_g_dabi = g_dabi("d,a,b,i").block(g_dabi_low, g_dabi_up);

                // block for t2_dcjk
                auto block_t2_dcjk = this_t2_dcjk("d,c,j,k");

                // block for g_cjkl
                auto block_g_cjkl = this_g_cjkl("c,j,k,l");

                // block for t2_abil
                block t2_abil_low{a_low, b_low, 0, 0};
                block t2_abil_up{a_up, b_up, n_tr_occ, n_tr_occ_inner};

                auto block_t2_abil = t2_right("a,b,i,l").block(t2_abil_low, t2_abil_up);

                t3("a,b,i,c,j,k") +=
                    (block_g_dabi * block_t2_dcjk) - (block_t2_abil * block_g_cjkl);
              }
            };

            // lambda function to compute v3
            auto compute_v3 = [&g_abij, &t1_this, n_tr_occ](std::size_t a,
                                                            std::size_t b,
                                                            std::size_t c, TArray &v3) {
              // index
              std::size_t a_low = a;
              std::size_t a_up = a + 1;
              std::size_t b_low = b;
              std::size_t b_up = b + 1;
              std::size_t c_low = c;
              std::size_t c_up = c + 1;

              {
                // block for g_bcjk
                block g_bcjk_low{b_low, c_low, 0, 0};
                block g_bcjk_up{b_up, c_up, n_tr_occ, n_tr_occ};

                auto block_g_bcjk = g_abij("b,c,j,k").block(g_bcjk_low, g_bcjk_up);

                // block for t1_ai
                block t1_ai_low{a_low, 0};
                block t1_ai_up{a_up, n_tr_occ};

                auto block_t_ai = t1_this("a,i").block(t1_ai_low, t1_ai_up);

                if (v3.is_initialized()) {
                  v3("b,c,j,k,a,i") += (block_g_bcjk * block_t_ai);
                } else {
                  v3("b,c,j,k,a,i") = (block_g_bcjk * block_t_ai);
                }
              }
            };

            // start (T) calculation
            double triple_energy = 0.0;

            // time spend in reduction
            double reduce_time = 0.0;
            // time spend in contraction
            double contraction_time = 0.0;
            // time spend in outer product
            double outer_product_time = 0.0;
            // time spend in permutation
            double permutation_time = 0.0;

            mpqc::time_point time00;
            mpqc::time_point time01;

            // local iteration
            std::size_t iter = 0;

            // global iteration
            std::size_t global_iter = 0;

            // make progress points
            const auto divide = 10;
            std::size_t increase = std::max(1.0, std::round(double(n_tr_vir) / divide));
            std::vector<std::size_t> progress_points;
            for (std::size_t i = 0; i < n_tr_vir; i += increase) {
              progress_points.push_back(i);
            }

            TA::set_default_world(this_world);

            // start loop over a, b, c
            for (auto a = 0; a < n_tr_vir; ++a) {
              std::size_t a_low = a;
              std::size_t a_up = a + 1;

              // block for g_ajkl
              block g_ajkl_low{a_low, 0, 0, 0};
              block g_ajkl_up{a_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};
              TArray g_ajkl;
              g_ajkl("a,j,k,l") = g_cjkl_global("a,j,k,l").block(g_ajkl_low, g_ajkl_up);

              // block for t2_dcjk
              block t2_dajk_low{0, a_low, 0, 0};
              block t2_dajk_up{n_tr_vir_inner, a_up, n_tr_occ, n_tr_occ};
              TArray t2_dajk;
              t2_dajk("d,a,j,k") = t2_left("d,a,j,k").block(t2_dajk_low, t2_dajk_up);

              for (auto b = 0; b <= a; ++b) {
                std::size_t b_low = b;
                std::size_t b_up = b + 1;

                // block for g_bjkl
                block g_bjkl_low{b_low, 0, 0, 0};
                block g_bjkl_up{b_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};
                TArray g_bjkl;
                g_bjkl("b,j,k,l") =
                    g_cjkl_global("b,j,k,l").block(g_bjkl_low, g_bjkl_up);

                // block for t2_dbjk
                block t2_dbjk_low{0, b_low, 0, 0};
                block t2_dbjk_up{n_tr_vir_inner, b_up, n_tr_occ, n_tr_occ};
                TArray t2_dbjk;
                t2_dbjk("d,b,j,k") = t2_left("d,b,j,k").block(t2_dbjk_low, t2_dbjk_up);

                for (auto c = 0; c <= b; ++c) {
                  global_iter++;

                  // round-robin distribute the loop
                  if (global_iter % size != rank) continue;

                  // inner loop
                  iter++;

                  // compute t3
                  TArray t3;
                  // abcijk contribution
                  // g^{da}_{bi}*t^{cd}_{kj} - g^{cj}_{kl}*t^{ab}_{il}

                  time00 = mpqc::now(this_world, accurate_time);
                  compute_t3_c(a, b, c, t3);
                  time01 = mpqc::now(this_world, accurate_time);
                  contraction_time += mpqc::duration_in_s(time00, time01);

                  // acbikj contribution
                  // g^{da}_{ci}*t^{db}_{kj} - g^{bk}_{jl}*t^{ac}_{il}
                  {
                    t3("a,c,i,b,k,j") = t3("a,b,i,c,j,k");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_t3(a, c, t3, t2_dbjk, g_bjkl);
                    time01 = mpqc::now(this_world, accurate_time);
                    contraction_time += mpqc::duration_in_s(time00, time01);
                  }

                  // cabkij contribution
                  // g^{dc}_{ak}*t^{db}_{ij} - g^{bi}_{jl}*t^{ca}_{kl}
                  {
                    t3("c,a,k,b,i,j") = t3("a,c,i,b,k,j");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_t3(c, a, t3, t2_dbjk, g_bjkl);
                    time01 = mpqc::now(this_world, accurate_time);
                    contraction_time += mpqc::duration_in_s(time00, time01);
                  }

                  // cbakji contribution
                  // g^{dc}_{bk}*t^{da}_{ji} - g^{aj}_{il}*t^{cb}_{kl}
                  {
                    t3("c,b,k,a,j,i") = t3("c,a,k,b,i,j");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_t3(c, b, t3, t2_dajk, g_ajkl);
                    time01 = mpqc::now(this_world, accurate_time);
                    contraction_time += mpqc::duration_in_s(time00, time01);
                  }

                  // bcajki contribution
                  // g^{db}_{cj}*t^{da}_{ki} - g^{ak}_{il}*t^{bc}_{jl}
                  {
                    t3("b,c,j,a,k,i") = t3("c,b,k,a,j,i");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_t3(b, c, t3, t2_dajk, g_ajkl);
                    time01 = mpqc::now(this_world, accurate_time);
                    contraction_time += mpqc::duration_in_s(time00, time01);
                  }

                  // bacjik contribution
                  // g^{db}_{aj}*t^{dc}_{ik} - g^{ci}_{kl}*t^{ba}_{jl}
                  {
                    t3("b,a,j,c,i,k") = t3("b,c,j,a,k,i");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_t3_c(b, a, c, t3);
                    time01 = mpqc::now(this_world, accurate_time);
                    contraction_time += mpqc::duration_in_s(time00, time01);
                  }

                  t3("a,b,c,i,j,k") = t3("b,a,j,c,i,k");
                  time00 = mpqc::now(this_world, accurate_time);
                  permutation_time += mpqc::duration_in_s(time01, time00);

                  // compute v3
                  TArray v3;

                  // bcajki contribution
                  // g^{bc}_{jk}*t^{a}_{i}
                  {
                    compute_v3(a, b, c, v3);
                    time01 = mpqc::now(this_world, accurate_time);
                    outer_product_time += mpqc::duration_in_s(time00, time01);
                  }

                  // acbikj contribution
                  // g^{ac}_{ik}*t^{b}_{j}
                  {
                    v3("a,c,i,k,b,j") = v3("b,c,j,k,a,i");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_v3(b, a, c, v3);
                    time01 = mpqc::now(this_world, accurate_time);
                    outer_product_time += mpqc::duration_in_s(time00, time01);
                  }

                  // abcijk contribution
                  // g^{ab}_{ij}*t^{c}_{k}
                  {
                    v3("a,b,i,j,c,k") = v3("a,c,i,k,b,j");
                    time00 = mpqc::now(this_world, accurate_time);
                    permutation_time += mpqc::duration_in_s(time01, time00);

                    compute_v3(c, a, b, v3);
                    time01 = mpqc::now(this_world, accurate_time);
                    outer_product_time += mpqc::duration_in_s(time00, time01);
                  }

                  v3("a,b,c,i,j,k") = v3("a,b,i,j,c,k");
                  time00 = mpqc::now(this_world, accurate_time);
                  permutation_time += mpqc::duration_in_s(time01, time00);

                  TArray result;
                  {
                    result("a,b,c,i,j,k") =
                        ((t3("a,b,c,i,j,k") + v3("a,b,c,i,j,k")) *
                            (4.0 * t3("a,b,c,i,j,k") + t3("a,b,c,k,i,j") +
                                t3("a,b,c,j,k,i") -
                                2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") +
                                    t3("a,b,c,j,i,k"))))
                            .set_world(this_world);
                    time01 = mpqc::now(this_world, accurate_time);
                    contraction_time += mpqc::duration_in_s(time00, time01);
                  }

                  // compute offset
                  std::size_t a_offset = tr_vir.tile(a).first;
                  std::size_t b_offset = tr_vir.tile(b).first;
                  std::size_t c_offset = tr_vir.tile(c).first;
                  std::array<std::size_t, 6> offset{
                      {a_offset, b_offset, c_offset, 0, 0, 0}};

                  double tmp_energy = 0.0;
                  if (b < a && c < b) {
                    time00 = mpqc::now(this_world, accurate_time);
                    auto ccsd_t_reduce = CCSD_T_Reduce(
                        this->orbital_energy(), trange1_engine->get_occ(),
                trange1_engine->get_nfrozen(), offset);
            tmp_energy = result("a,b,c,i,j,k").reduce(ccsd_t_reduce);
            time01 = mpqc::now(this_world, accurate_time);
            reduce_time += mpqc::duration_in_s(time00, time01);
            tmp_energy *= 2;
          }
          // boundary condition
          else {
            time00 = mpqc::now(this_world, accurate_time);
            auto ccsd_t_reduce = CCSD_T_ReduceSymm(
                this->orbital_energy(), trange1_engine->get_occ(),
                trange1_engine->get_nfrozen(), offset);
            tmp_energy = result("a,b,c,i,j,k").reduce(ccsd_t_reduce);
            time01 = mpqc::now(this_world, accurate_time);
            reduce_time += mpqc::duration_in_s(time00, time01);
          }

          triple_energy += tmp_energy;
        }  // loop of c
      }    // loop of b

      // print the progress
      if (rank == 0) {
        util::print_progress(a, a + 1, progress_points);
      }
    }  // loop of a
    this_world.gop.fence();
    global_world.gop.fence();

    TA::set_default_world(global_world);

    if (this->verbose()) {
      // loop over all rank and print
      for (auto i = 0; i < size; ++i) {
        global_world.gop.fence();
        if (rank == i) {
          // print out process n
          std::cout << "Process " << rank << " Time: " << std::endl;
          std::cout << "Iter: " << iter << std::endl;
          std::cout << "Permutation Time: " << permutation_time << " S"
                    << std::endl;
          std::cout << "Contraction Time: " << contraction_time << " S"
                    << std::endl;
          std::cout << "Outer Product Time: " << outer_product_time << " S"
                    << std::endl;
          std::cout << "Reduce Time: " << reduce_time << " S" << std::endl
                    << std::endl;
        }
      }
    }

    // print out all process time
    global_world.gop.sum(iter);
    global_world.gop.sum(permutation_time);
    global_world.gop.sum(contraction_time);
    global_world.gop.sum(outer_product_time);
    global_world.gop.sum(reduce_time);

    ExEnv::out0() << "Process All Time: " << std::endl;
    ExEnv::out0() << "Iter: " << iter << std::endl;
    ExEnv::out0() << "Permutation Time: " << permutation_time << " S"
                  << std::endl;
    ExEnv::out0() << "Contraction Time: " << contraction_time << " S"
                  << std::endl;
    ExEnv::out0() << "Outer Product Time: " << outer_product_time << " S"
                  << std::endl;
    ExEnv::out0() << "Reduce Time: " << reduce_time << " S" << std::endl
                  << std::endl;

    global_world.gop.sum(triple_energy);

    //    ExEnv::out0() << "(T) Energy: " << triple_energy << std::endl
    //                  << std::endl;

    // manually clean replicated array
    if (size > 1) {
      t1_this = TArray();
      g_cjkl_global = TArray();

      // warning! this line is important, it make sure all Array object was
      // cleaned before destruction of this_world
      TArray::wait_for_lazy_cleanup(this_world);
      global_world.gop.fence();
      world_ptr.reset();
    }

    global_world.gop.fence();
    return triple_energy;
  }

  double compute_ccsd_t_fine_grain(TArray &t1, TArray &t2) {
    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    // get integral
    TArray g_cjkl = get_aijk();
    TArray g_abij = get_abij();
    TArray g_dabi = get_abci();

    // T2
    TArray t2_left = t2;
    TArray t2_right = t2;

    if (reblock_inner_) {
      reblock_inner_t2(t2_left, t2_right);
    }

    auto trange1_engine = this->local_trange1_engine() == nullptr
                              ? this->trange1_engine()
                              : this->local_trange1_engine();

    // get trange1
    auto tr_occ = trange1_engine->get_active_occ_tr1();
    auto tr_vir = trange1_engine->get_vir_tr1();

    auto n_tr_occ = trange1_engine->get_active_occ_blocks();
    auto n_tr_vir = trange1_engine->get_vir_blocks();
    auto n_tr_occ_inner = n_tr_occ;
    auto n_tr_vir_inner = n_tr_vir;
    if (reblock_inner_) {
      n_tr_occ_inner = tr_occ_inner_.tiles_range().second;
      n_tr_vir_inner = tr_vir_inner_.tiles_range().second;
    }

    // lambda function to compute t3
    auto compute_t3 = [&](std::size_t a_low, std::size_t b_low,
                          std::size_t c_low, std::size_t a_up, std::size_t b_up,
                          std::size_t c_up, TArray &t3) {

      typedef std::vector<std::size_t> block;

      {
        // block for g_dbai
        block g_dabi_low{0, a_low, b_low, 0};
        block g_dabi_up{n_tr_vir_inner, a_up, b_up, n_tr_occ};

        auto block_g_dabi = g_dabi("d,a,b,i").block(g_dabi_low, g_dabi_up);

        // block for t2_cdk
        block t2_dcjk_low{0, c_low, 0, 0};
        block t2_dcjk_up{n_tr_vir_inner, c_up, n_tr_occ, n_tr_occ};

        auto block_t2_dcjk = t2_left("d,c,j,k").block(t2_dcjk_low, t2_dcjk_up);

        // block for g_cjkl
        block g_cjkl_low{c_low, 0, 0, 0};
        block g_cjkl_up{c_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

        auto block_g_cjkl = g_cjkl("c,j,k,l").block(g_cjkl_low, g_cjkl_up);

        // block for t2_abil
        block t2_abil_low{a_low, b_low, 0, 0};
        block t2_abil_up{a_up, b_up, n_tr_occ, n_tr_occ_inner};

        auto block_t2_abil = t2_right("a,b,i,l").block(t2_abil_low, t2_abil_up);

        if (t3.is_initialized()) {
          t3("a,b,i,c,j,k") +=
              block_g_dabi * block_t2_dcjk - block_t2_abil * block_g_cjkl;
        } else {
          t3("a,b,i,c,j,k") =
              block_g_dabi * block_t2_dcjk - block_t2_abil * block_g_cjkl;
        }
      }
    };

    // lambda function to compute v3
    auto compute_v3 = [&](std::size_t a_low, std::size_t b_low,
                          std::size_t c_low, std::size_t a_up, std::size_t b_up,
                          std::size_t c_up, TArray &v3) {
      typedef std::vector<std::size_t> block;

      {
        // block for g_bcjk
        block g_bcjk_low{b_low, c_low, 0, 0};
        block g_bcjk_up{b_up, c_up, n_tr_occ, n_tr_occ};

        auto block_g_bcjk = g_abij("b,c,j,k").block(g_bcjk_low, g_bcjk_up);

        // block for t1_ai
        block t1_ai_low{a_low, 0};
        block t1_ai_up{a_up, n_tr_occ};

        auto block_t_ai = t1("a,i").block(t1_ai_low, t1_ai_up);

        if (v3.is_initialized()) {
          v3("b,c,j,k,a,i") += block_g_bcjk * block_t_ai;
        } else {
          v3("b,c,j,k,a,i") = block_g_bcjk * block_t_ai;
        }
      }
    };

    double triple_energy = 0.0;

    // number of unocc blocks to load at each iteration
    std::size_t increase = increase_;
    if (increase > n_tr_vir) {
      increase = n_tr_vir;
    }
    std::size_t a_increase = increase;
    std::size_t b_increase = increase;
    std::size_t c_increase = increase;

    std::size_t occ_block_size =
        trange1_engine->get_occ_block_size();
    std::size_t vir_block_size =
        trange1_engine->get_vir_block_size();
    std::size_t n_blocks =
        increase * increase * increase * n_tr_occ * n_tr_occ * n_tr_occ;
    double mem = (n_blocks * std::pow(occ_block_size, 3) *
                  std::pow(vir_block_size, 3) * 8) /
                 1.0e9;

    if (t1.world().rank() == 0) {
      std::cout << "Increase in the loop " << increase << std::endl;
      std::cout << "Number of blocks at each iteration " << n_blocks
                << std::endl;
      std::cout << "Size of T3 or V3 at each iteration " << mem << " GB"
                << std::endl;
    }

    double t3_time = 0.0;
    double v3_time = 0.0;
    double reduce_time = 0.0;
    double t3_permute_time = 0.0;
    double v3_permute_time = 0.0;
    double contraction_time1 = 0.0;
    double contraction_time2 = 0.0;
    double contraction_time3 = 0.0;
    double contraction_time4 = 0.0;
    double contraction_time5 = 0.0;
    double contraction_time6 = 0.0;
    double v3_contraction_time = 0.0;
    mpqc::time_point time0;
    mpqc::time_point time1;
    mpqc::time_point time2;
    mpqc::time_point time3;
    mpqc::time_point time00;
    mpqc::time_point time01;

    // index in virtual blocks
    std::size_t a = 0;
    std::size_t b = 0;
    std::size_t c = 0;

    // number of blocks computed
    std::size_t n_blocks_computed = 0;

    const auto divide = 10;
    std::size_t progress_increase =
        std::max(1.0, std::round(double(n_tr_vir) / divide));
    std::vector<std::size_t> progress_points;
    for (int i = 0; i < n_tr_vir; i += progress_increase) {
      progress_points.push_back(i);
    }
    // loop over virtual blocks
    while (a < n_tr_vir) {
      b = 0;
      if (a + increase >= n_tr_vir) {
        a_increase = (n_tr_vir - a);
      } else {
        a_increase = increase;
      }

      std::size_t a_end = a + a_increase - 1;
      while (b <= a_end) {
        c = 0;

        if (b + increase - 1 > a_end) {
          if (b == a_end) {
            b_increase = 1;
          } else {
            b_increase = (a_end - b);
          }
        } else {
          b_increase = increase;
        }

        std::size_t b_end = b + b_increase - 1;

        while (c <= b_end) {
          if (c + increase - 1 > b_end) {
            if (c == b_end) {
              c_increase = 1;
            } else {
              c_increase = (b_end - c);
            }
          } else {
            c_increase = increase;
          }

          std::size_t c_end = c + c_increase - 1;

          //                            std::cout << a << " " << b << " " << c
          //                            << std::endl;
          std::size_t a_low = a;
          std::size_t a_up = a + a_increase;
          std::size_t b_low = b;
          std::size_t b_up = b + b_increase;
          std::size_t c_low = c;
          std::size_t c_up = c + c_increase;

          std::size_t blocks = (a_up - a_low) * (b_up - b_low) *
                               (c_up - c_low) * n_tr_occ * n_tr_occ * n_tr_occ;

          n_blocks_computed += blocks;

          time0 = mpqc::now(world, accurate_time);
          // compute t3
          TArray t3;
          // abcijk contribution
          // g^{da}_{bi}*t^{cd}_{kj} - g^{cj}_{kl}*t^{ab}_{il}
          time00 = mpqc::now(world, accurate_time);
          compute_t3(a_low, b_low, c_low, a_up, b_up, c_up, t3);
          time01 = mpqc::now(world, accurate_time);
          contraction_time1 += mpqc::duration_in_s(time00, time01);

          // acbikj contribution
          // g^{da}_{ci}*t^{db}_{kj} - g^{bk}_{jl}*t^{ac}_{il}
          {
            t3("a,c,i,b,k,j") = t3("a,b,i,c,j,k");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_t3(a_low, c_low, b_low, a_up, c_up, b_up, t3);
            time01 = mpqc::now(world, accurate_time);
            contraction_time2 += mpqc::duration_in_s(time00, time01);
          }

          // cabkij contribution
          // g^{dc}_{ak}*t^{db}_{ij} - g^{bi}_{jl}*t^{ca}_{kl}
          {
            t3("c,a,k,b,i,j") = t3("a,c,i,b,k,j");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_t3(c_low, a_low, b_low, c_up, a_up, b_up, t3);
            time01 = mpqc::now(world, accurate_time);
            contraction_time3 += mpqc::duration_in_s(time00, time01);
          }

          // cbakji contribution
          // g^{dc}_{bk}*t^{da}_{ji} - g^{aj}_{il}*t^{cb}_{kl}
          {
            t3("c,b,k,a,j,i") = t3("c,a,k,b,i,j");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_t3(c_low, b_low, a_low, c_up, b_up, a_up, t3);
            time01 = mpqc::now(world, accurate_time);
            contraction_time4 += mpqc::duration_in_s(time00, time01);
          }

          // bcajki contribution
          // g^{db}_{cj}*t^{da}_{ki} - g^{ak}_{il}*t^{bc}_{jl}
          {
            t3("b,c,j,a,k,i") = t3("c,b,k,a,j,i");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_t3(b_low, c_low, a_low, b_up, c_up, a_up, t3);
            time01 = mpqc::now(world, accurate_time);
            contraction_time5 += mpqc::duration_in_s(time00, time01);
          }

          // bacjik contribution
          // g^{db}_{aj}*t^{dc}_{ik} - g^{ci}_{kl}*t^{ba}_{jl}
          {
            t3("b,a,j,c,i,k") = t3("b,c,j,a,k,i");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_t3(b_low, a_low, c_low, b_up, a_up, c_up, t3);
            time01 = mpqc::now(world, accurate_time);
            contraction_time6 += mpqc::duration_in_s(time00, time01);
          }

          time00 = mpqc::now(world, accurate_time);
          t3("a,b,c,i,j,k") = t3("b,a,j,c,i,k");
          time01 = mpqc::now(world, accurate_time);
          t3_permute_time += mpqc::duration_in_s(time00, time00);

          time1 = mpqc::now(world, accurate_time);
          t3_time += mpqc::duration_in_s(time0, time1);

          // compute v3
          TArray v3;

          // bcajki contribution
          // g^{bc}_{jk}*t^{a}_{i}
          {
            time00 = mpqc::now(world, accurate_time);
            compute_v3(a_low, b_low, c_low, a_up, b_up, c_up, v3);
            time01 = mpqc::now(world, accurate_time);
            v3_contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // acbikj contribution
          // g^{ac}_{ik}*t^{b}_{j}
          {
            v3("a,c,i,k,b,j") = v3("b,c,j,k,a,i");
            time00 = mpqc::now(world, accurate_time);
            v3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_v3(b_low, a_low, c_low, b_up, a_up, c_up, v3);
            time01 = mpqc::now(world, accurate_time);
            v3_contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // abcijk contribution
          // g^{ab}_{ij}*t^{c}_{k}
          {
            v3("a,b,i,j,c,k") = v3("a,c,i,k,b,j");
            time00 = mpqc::now(world, accurate_time);
            v3_permute_time += mpqc::duration_in_s(time01, time00);

            compute_v3(c_low, a_low, b_low, c_up, a_up, b_up, v3);
            time01 = mpqc::now(world, accurate_time);
            v3_contraction_time += mpqc::duration_in_s(time00, time01);
          }

          time00 = mpqc::now(world, accurate_time);
          v3("a,b,c,i,j,k") = v3("a,b,i,j,c,k");
          time01 = mpqc::now(world, accurate_time);
          v3_permute_time += mpqc::duration_in_s(time00, time01);

          time2 = mpqc::now(world, accurate_time);
          v3_time += mpqc::duration_in_s(time1, time2);

          // compute offset
          std::size_t a_offset = tr_vir.tile(a).first;
          std::size_t b_offset = tr_vir.tile(b).first;
          std::size_t c_offset = tr_vir.tile(c).first;
          //                            std::cout << a_offset << " " <<
          //                            b_offset
          //                            << " " << c_offset << std::endl;
          std::array<std::size_t, 6> offset{
              {a_offset, b_offset, c_offset, 0, 0, 0}};

          double tmp_energy = 0.0;
          if (b_end < a && c_end < b) {
            auto ccsd_t_reduce = CCSD_T_Reduce(
                this->orbital_energy(), trange1_engine->get_occ(),
                trange1_engine->get_nfrozen(), offset);
            tmp_energy = ((t3("a,b,c,i,j,k") + v3("a,b,c,i,j,k")) *
                          (4.0 * t3("a,b,c,i,j,k") + t3("a,b,c,k,i,j") +
                           t3("a,b,c,j,k,i") -
                           2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") +
                                t3("a,b,c,j,i,k"))))
                             .reduce(ccsd_t_reduce);

            tmp_energy *= 2;
          } else {
            auto ccsd_t_reduce = CCSD_T_ReduceSymm(
                this->orbital_energy(), trange1_engine->get_occ(),
                trange1_engine->get_nfrozen(), offset);
            tmp_energy = ((t3("a,b,c,i,j,k") + v3("a,b,c,i,j,k")) *
                          (4.0 * t3("a,b,c,i,j,k") + t3("a,b,c,k,i,j") +
                           t3("a,b,c,j,k,i") -
                           2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") +
                                t3("a,b,c,j,i,k"))))
                             .reduce(ccsd_t_reduce);
          }

          time3 = mpqc::now(world, accurate_time);
          reduce_time += mpqc::duration_in_s(time2, time3);

          triple_energy += tmp_energy;

          c += c_increase;
        }
        b += b_increase;
      }

      if (t1.world().rank() == 0) {
        util::print_progress(a, a + increase, progress_points);
      }
      a += a_increase;
    }

    if (t1.world().rank() == 0) {
      std::cout << "Total Blocks Computed  " << n_blocks_computed;
      std::cout << " from " << std::pow(n_tr_occ, 3) * std::pow(n_tr_vir, 3)
                << std::endl;

      std::cout << "T3 Total Time: " << t3_time << " S \n";
      std::cout << "T3 Permutation Time: " << t3_permute_time << " S \n";
      std::cout << "T3 Contraction Time1: " << contraction_time1 << " S \n";
      std::cout << "T3 Contraction Time2: " << contraction_time2 << " S \n";
      std::cout << "T3 Contraction Time3: " << contraction_time3 << " S \n";
      std::cout << "T3 Contraction Time4: " << contraction_time4 << " S \n";
      std::cout << "T3 Contraction Time5: " << contraction_time5 << " S \n";
      std::cout << "T3 Contraction Time6: " << contraction_time6 << " S \n";
      std::cout << "V3 Total Time: " << v3_time << " S \n";
      std::cout << "V3 Contraction Time: " << v3_contraction_time << " S \n";
      std::cout << "V3 Permutation Time: " << v3_permute_time << " S \n";
      std::cout << "Reduction Total Time: " << reduce_time << " S \n";
    }
    return triple_energy;
  }

  // compute and store all t3 amplitudes, not recommended for performance
  // computing
  double compute_ccsd_t_straight(const TArray &t1, const TArray &t2) {
    // get integral
    TArray g_cjkl = get_aijk();
    TArray g_dabi = get_abci();
    TArray g_abij = get_abij();

    TArray t2_left = t2;
    TArray t2_right = t2;

    if (reblock_inner_) {
      reblock_inner_t2(t2_left, t2_right);
    }

    // compute t3
    TArray t3;
    t3("a,b,c,i,j,k") = g_dabi("d,a,b,i") * t2_left("d,c,j,k") -
                        g_cjkl("c,j,k,l") * t2_right("a,b,i,l");
    t3("a,b,c,i,j,k") = t3("a,b,c,i,j,k") + t3("a,c,b,i,k,j") +
                        t3("c,a,b,k,i,j") + t3("c,b,a,k,j,i") +
                        t3("b,c,a,j,k,i") + t3("b,a,c,j,i,k");

    // compute v3
    TArray v3;
    v3("a,b,c,i,j,k") = g_abij("a,b,i,j") * t1("c,k");
    v3("a,b,c,i,j,k") =
        v3("a,b,c,i,j,k") + v3("b,c,a,j,k,i") + v3("a,c,b,i,k,j");

    std::array<std::size_t, 6> offset{{0, 0, 0, 0, 0, 0}};
    auto trange1_engine = this->local_trange1_engine() == nullptr
                          ? this->trange1_engine()
                          : this->local_trange1_engine();

    auto ccsd_t_reduce = CCSD_T_Reduce(
        this->orbital_energy(), trange1_engine->get_occ(),
        trange1_engine->get_nfrozen(), offset);

    double triple_energy =
        ((t3("a,b,c,i,j,k") + v3("a,b,c,i,j,k")) *
         (4.0 * t3("a,b,c,i,j,k") + t3("a,b,c,k,i,j") + t3("a,b,c,j,k,i") -
          2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") + t3("a,b,c,j,i,k"))))
            .reduce(ccsd_t_reduce);
    triple_energy = triple_energy / 3.0;
    return triple_energy;
  }

  // performs Laplace transform perturbative triple correction to CCSD energy
  double compute_ccsd_t_laplace_transform(const TArray &t1, const TArray &t2) {
    // get integral
    TArray g_cjkl = get_aijk();
    TArray g_dabi = get_abci();
    TArray g_abij = get_abij();

    // obtaining DF-integrals
    TArray Xab;
    TArray Xai;
    if (this->is_df()) {
      /// get three center integral (X|ab)
      Xab = this->lcao_factory().compute(L"(Κ|G|a b)[inv_sqr]");

      // TArray Xai;
      Xai = this->lcao_factory().compute(L"(Κ|G|a i)[inv_sqr]");
    }

    // definition of orbital spaces
    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    // copying orbital energies into Eigen vector
    Eigen::VectorXd const &e_orb = *this->orbital_energy();

    // definition of alpha parameter required for the Gaussian quadrature
    // defined as in paper: Constans, Pere, Philippe Y. Ayala, and Gustavo E.
    // Scuseria.
    //"Scaling reduction of the perturbative triples correction (T) to coupled
    // cluster
    // theory via Laplace transform formalism." The Journal of Chemical Physics
    // 113.23 (2000): 10451-10458.
    double alpha = 3.0 * (e_orb(n_occ) - e_orb(n_occ - 1));

    // definition of number of quadrature points
    // should be included in input as a parameter
    int n = n_laplace_quad_;

    // defining the weights and roots for quadrature
    Eigen::VectorXd x(n);
    Eigen::VectorXd w(n);
    // function that evaluates Gaussian weights and roots for quadrature using
    // orthogonal polynomials
    mpqc::math::gauss_legendre(n, w, x);

    double triple_energy = 0.0;

    // get trange1
    auto n_tr_occ = this->trange1_engine()->get_active_occ_blocks();
    auto n_tr_vir = this->trange1_engine()->get_vir_blocks();

    auto &global_world = this->wfn_world()->world();
    // split global_world
    const auto rank = global_world.rank();
    const auto size = global_world.size();

    // required for coarse grain parallelism
    madness::World *tmp_ptr;
    std::shared_ptr<madness::World> world_ptr;

    if (size > 1) {
      SafeMPI::Group group = global_world.mpi.comm().Get_group().Incl(1, &rank);
      SafeMPI::Intracomm comm = global_world.mpi.comm().Create(group);
      world_ptr = std::make_shared<madness::World>(comm);
      tmp_ptr = world_ptr.get();
    } else {
      tmp_ptr = &global_world;
    }
    auto &this_world = *tmp_ptr;
    global_world.gop.fence();

    // loop over number of quadrature points
    for (auto m = 0; m < n; m++) {
      // conversion of integrals  and amplitudes to Laplace-transform form.
      TArray g_dabi_lt = g_dabi_laplace_transform(
          g_dabi, *this->orbital_energy(), n_occ, n_frozen, x(m));
      TArray t2_oou_lt = t2_oou_laplace_transform(t2, *this->orbital_energy(),
                                                  n_occ, n_frozen, x(m));

      TArray g_cjkl_lt = g_cjkl_laplace_transform(
          g_cjkl, *this->orbital_energy(), n_occ, n_frozen, x(m));
      TArray t2_ouu_lt = t2_ouu_laplace_transform(t2, *this->orbital_energy(),
                                                  n_occ, n_frozen, x(m));

      // conversion of the 3-center integrals required in DF into
      // Laplace-transform form
      TArray Xab_lt;
      TArray Xai_lt;
      if (this->is_df()) {
        Xab_lt = Xab_laplace_transform(Xab, *this->orbital_energy(), n_occ,
                                       n_frozen, x(m));

        Xai_lt = Xai_laplace_transform(Xai, *this->orbital_energy(), n_occ,
                                       n_frozen, x(m));
      }

      this->wfn_world()->world().gop.fence();

      // energy contribution for each quadrature point of various intermediates
      // as defined in "Scaling reduction of the perturbative triples correction
      // (T) to coupled cluster theory via Laplace transform formalism"
      // energy_vv arise from unoccupied-unoccupied contribution obtained from
      // W*W from Eq. (5). (Loop over "f")
      // energy_ov arise from occupied-occupied contribution obtained from W*W
      // from Eq. (5). (Loop over "m")
      // energy_vo arise from unoccupied-occupied contribution obtained from W*W
      // from Eq. (5). (Loop over "f" and "m")
      // energy_v_s arise from unoccupied-V contribution (W*V). (Loop over "f")
      // energy_o_s arise from occupied-V contribution (W*V). (Loop over "m")
      double energy_vv = 0.0;
      double energy_oo = 0.0;
      double energy_vo = 0.0;
      double energy_v_s = 0.0;
      double energy_o_s = 0.0;

      // intermediate required for OV5 computation within DF formalism
      TArray gg1;
      if (this->is_df()) {
        gg1("X,Y") = Xai_lt("X,c,i") * Xai_lt("Y,c,i");
      }

      double E_OV5 = 0;
      double E_O2V4_vo = 0;

      // blocking over a & b to avoid all intermediate arrays over 3 and 4
      // unoccupied indices
      TA::set_default_world(this_world);
      std::size_t global_iter = 0;
      for (auto a = 0; a < n_tr_vir; ++a) {
        for (auto b = 0; b < n_tr_vir; ++b) {
          ++global_iter;

          if (global_iter % size != rank) continue;

          std::size_t a_low = a;
          std::size_t a_up = a + 1;
          std::size_t b_low = b;
          std::size_t b_up = b + 1;

          typedef std::vector<std::size_t> block;

          // blocking g_dabi over b
          block g_dabi_low_b{0, b_low, 0, 0};
          block g_dabi_up_b{n_tr_vir, b_up, n_tr_vir, n_tr_occ};

          auto block_g_dabi_lt_b =
              g_dabi_lt("e,b,c,i").block(g_dabi_low_b, g_dabi_up_b);

          // blocking g_dabi over a
          block g_dabi_low_a{0, a_low, 0, 0};
          block g_dabi_up_a{n_tr_vir, a_up, n_tr_vir, n_tr_occ};
          auto block_g_dabi_lt_a =
              g_dabi_lt("f,a,c,i").block(g_dabi_low_a, g_dabi_up_a);

          // blocking t2 over a (required for both T1 and T2 intermediates)
          block t2_oou_lt_low_a{0, a_low, 0, 0};
          block t2_oou_lt_up_a{n_tr_vir, a_up, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_a_T12 =
              t2_oou_lt("e,a,i,j").block(t2_oou_lt_low_a, t2_oou_lt_up_a);

          // blocking t2 over b (required for T2 intermediate)
          block t2_oou_lt_low_b{0, b_low, 0, 0};
          block t2_oou_lt_up_b{n_tr_vir, b_up, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_b_T2 =
              t2_oou_lt("f,b,j,i").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

          // blocking t2 over b (required for T1 intermediate)
          auto block_t2_oou_lt_b_T1 =
              t2_oou_lt("f,b,i,j").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

          TArray T1;
          TArray T2;
          // intermediate "amplitudes" required for the OV5 step.
          // Computational cost of T1 & T2 is O2V4
          T1("e,b,a,f") = block_t2_oou_lt_a_T12 * block_t2_oou_lt_b_T1;
          T2("e,b,a,f") = block_t2_oou_lt_a_T12 * block_t2_oou_lt_b_T2;

          {
            // Computational cost of OV5
            TArray G;
            G("e,b,a,f") = block_g_dabi_lt_b * block_g_dabi_lt_a;

            E_OV5 +=
                TA::dot((G("e,b,a,f")), (T2("e,b,a,f") - 2.0 * T1("e,b,a,f")));
          }
          if (this->is_df()) {
            // DF approximation reduces computational cost OV5 to at most AV4
            // (A=Aux basis set)
            auto n_tr_aux = Xab.range().upbound()[0];

            block Xab_low_e{0, 0, b_low};
            block Xab_up_e{n_tr_aux, n_tr_vir, b_up};
            auto block_Xab_lt_e = Xab_lt("X,e,b").block(Xab_low_e, Xab_up_e);

            block Xab_low_f{0, 0, a_low};
            block Xab_up_f{n_tr_aux, n_tr_vir, a_up};
            auto block_Xab_lt_f = Xab_lt("Y,f,a").block(Xab_low_f, Xab_up_f);

            TArray gT1;
            gT1("a,f,X") =
                block_Xab_lt_e * (4.0 * T2("e,b,a,f") - 2.0 * T1("e,b,a,f"));

            TArray ggT1;
            ggT1("X,Y") = gT1("a,f,X") * block_Xab_lt_f;
            E_OV5 += TA::dot((ggT1("X,Y")), (gg1("X,Y")));

            TArray gg2;
            gg2("e,b,Y") = block_g_dabi_lt_b * Xai_lt("Y,c,i");

            TArray ggT2;
            ggT2("a,f,Y") =
                gg2("e,b,Y") * (4.0 * T2("e,b,a,f") - 2.0 * T1("e,b,a,f"));

            E_OV5 -= TA::dot((ggT2("a,f,Y")), (block_Xab_lt_f));
          } else {
            TArray G;

            block g_dabi_low_ca{0, 0, a_low, 0};
            block g_dabi_up_ca{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
            auto block_g_dabi_lt_ca =
                g_dabi_lt("f,c,a,i").block(g_dabi_low_ca, g_dabi_up_ca);

            block g_dabi_low_cb{0, 0, b_low, 0};
            block g_dabi_up_cb{n_tr_vir, n_tr_vir, b_up, n_tr_occ};
            auto block_g_dabi_lt_cb =
                g_dabi_lt("e,c,b,i").block(g_dabi_low_cb, g_dabi_up_cb);

            G("e,b,a,f") = block_g_dabi_lt_cb * block_g_dabi_lt_ca -
                           block_g_dabi_lt_b * block_g_dabi_lt_ca;

            E_OV5 += TA::dot((G("e,b,a,f")),
                             (4.0 * T2("e,b,a,f") - 2.0 * T1("e,b,a,f")));
          }

          // Mixed term contributions
          {
            TArray T1;
            TArray T2;

            auto block_t2_oou_lt_eb =
                t2_oou_lt("e,b,i,j").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

            block g_cjkl_low_a{a_low, 0, 0, 0};
            block g_cjkl_up_a{a_up, n_tr_occ, n_tr_occ, n_tr_occ};

            auto block_g_dabi_lt_bji =
                g_cjkl_lt("a,j,i,n").block(g_cjkl_low_a, g_cjkl_up_a);
            auto block_g_dabi_lt_bij =
                g_cjkl_lt("a,i,j,n").block(g_cjkl_low_a, g_cjkl_up_a);

            T1("e,a,b,n") = block_t2_oou_lt_eb * block_g_dabi_lt_bji;
            T2("e,a,b,n") = block_t2_oou_lt_eb * block_g_dabi_lt_bij;

            if (this->is_df()) {
              auto n_tr_aux = Xab.range().upbound()[0];

              TArray G1;
              TArray G2;
              TArray G3;

              auto block_t2_ouu_lt_cb =
                  t2_ouu_lt("c,b,i,n").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

              block t2_oou_lt_low_bc{b_low, 0, 0, 0};
              block t2_oou_lt_up_bc{b_up, n_tr_vir, n_tr_occ, n_tr_occ};
              auto block_t2_ouu_lt_bc =
                  t2_ouu_lt("b,c,i,n").block(t2_oou_lt_low_bc, t2_oou_lt_up_bc);

              auto block_g_dabi_lt_ac =
                  g_dabi_lt("e,a,c,i").block(g_dabi_low_a, g_dabi_up_a);

              block g_dabi_low_ca{0, 0, a_low, 0};
              block g_dabi_up_ca{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
              auto block_g_dabi_lt_ca =
                  g_dabi_lt("e,c,a,i").block(g_dabi_low_ca, g_dabi_up_ca);

              G2("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_cb;
              E_O2V4_vo += TA::dot((-G2("e,c,a,i")),
                                   (2.0 * T1("e,c,a,i") - T2("e,c,a,i")));

              TArray G;
              G("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_bc;
              E_O2V4_vo +=
                  TA::dot(G("e,c,a,i"), (-2.0 * T2("e,c,a,i") + T1("e,c,a,i")));

              TArray gt;
              gt("b,n,X") = Xai_lt("X,c,i") *
                            (2.0 * block_t2_ouu_lt_cb - block_t2_ouu_lt_bc);

              TArray gtT;
              gtT("e,a,X") =
                  gt("b,n,X") * (2.0 * T1("e,a,b,n") - T2("e,a,b,n"));

              block Xab_low{0, 0, a_low};
              block Xab_up{n_tr_aux, n_tr_vir, a_up};
              auto block_Xab_lt = Xab_lt("X,e,a").block(Xab_low, Xab_up);
              E_O2V4_vo += TA::dot(gtT("e,a,X"), block_Xab_lt);

            } else {
              TArray G1;
              TArray G2;
              TArray G3;

              auto block_t2_ouu_lt_cb =
                  t2_ouu_lt("c,b,i,n").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

              block t2_oou_lt_low_bc{b_low, 0, 0, 0};
              block t2_oou_lt_up_bc{b_up, n_tr_vir, n_tr_occ, n_tr_occ};
              auto block_t2_ouu_lt_bc =
                  t2_ouu_lt("b,c,i,n").block(t2_oou_lt_low_bc, t2_oou_lt_up_bc);

              auto block_g_dabi_lt_ac =
                  g_dabi_lt("e,a,c,i").block(g_dabi_low_a, g_dabi_up_a);

              block g_dabi_low_ca{0, 0, a_low, 0};
              block g_dabi_up_ca{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
              auto block_g_dabi_lt_ca =
                  g_dabi_lt("e,c,a,i").block(g_dabi_low_ca, g_dabi_up_ca);

              G1("e,a,b,n") = block_g_dabi_lt_ca * block_t2_ouu_lt_cb;
              G2("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_cb;
              G3("e,a,b,n") = block_g_dabi_lt_ca * block_t2_ouu_lt_bc;
              E_O2V4_vo +=
                  TA::dot(2.0 * G1("e,c,a,i") - G2("e,c,a,i") - G3("e,c,a,i"),
                          (2.0 * T1("e,c,a,i") - T2("e,c,a,i")));

              TArray G;
              G("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_bc;
              E_O2V4_vo +=
                  TA::dot(G("e,c,a,i"), (-2.0 * T2("e,c,a,i") + T1("e,c,a,i")));
            }
          }
        }
      }
      TA::set_default_world(global_world);
      global_world.gop.sum(E_OV5);
      global_world.gop.sum(E_O2V4_vo);

      // intermediates required for O2V4 computation within DF formalism
      TArray gt1;
      TArray gt2;
      if (this->is_df()) {
        gt1("f,X,i") = Xai_lt("X,b,j") * t2_oou_lt("f,b,j,i");
        gt2("f,X,i") = Xai_lt("X,b,j") * t2_oou_lt("f,b,i,j");
      }

      // blocking over a to avoid all intermediate arrays over 3 unoccupied
      // indices
      double E_O2V4_vv = 0;
      TA::set_default_world(this_world);
      global_iter = 0;
      for (auto a = 0; a < n_tr_vir; ++a) {
        ++global_iter;
        if (global_iter % size != rank) continue;

        std::size_t a_low = a;
        std::size_t a_up = a + 1;

        typedef std::vector<std::size_t> block;
        {
          // blocking g_dabi over a & b
          block g_dabi_low_ab{0, a_low, 0, 0};
          block g_dabi_up_ab{n_tr_vir, a_up, n_tr_vir, n_tr_occ};
          auto block_g_dabi_lt_ab =
              g_dabi_lt("e,a,b,j").block(g_dabi_low_ab, g_dabi_up_ab);

          TArray G2;
          G2("e,a,i,f") = block_g_dabi_lt_ab * t2_oou_lt("f,b,j,i");

          TArray G3;
          G3("e,a,i,f") = block_g_dabi_lt_ab * t2_oou_lt("f,b,i,j");

          E_O2V4_vv +=
              TA::dot(G2("e,a,i,f"), G2("f,a,i,e") - 4.0 * G3("f,a,i,e"));
          {
            if (this->is_df()) {
              auto n_tr_aux = Xab.range().upbound()[0];

              block Xab_low{0, 0, a_low};
              block Xab_up{n_tr_aux, n_tr_vir, a_up};
              auto block_Xab_lt_ba = Xab_lt("X,e,a").block(Xab_low, Xab_up);

              TArray G4;
              G4("e,a,i,f") = block_Xab_lt_ba * gt1("f,X,i");
              E_O2V4_vv +=
                  TA::dot(G4("e,a,i,f"), G4("f,a,i,e") - 4.0 * G2("f,a,i,e"));
              E_O2V4_vv +=
                  TA::dot(G3("e,a,i,f"), 2.0 * G4("f,a,i,e") + G3("f,a,i,e"));

              TArray G1;
              G1("e,a,i,f") = block_Xab_lt_ba * gt2("f,X,i");

              E_O2V4_vv += TA::dot(
                  G1("e,a,i,f"), 8.0 * G2("f,a,i,e") + 4.0 * G1("f,a,i,e") -
                                     4.0 * G4("f,a,i,e") - 4.0 * G3("f,a,i,e"));

            } else {
              block g_dabi_low_ba{0, 0, a_low, 0};
              block g_dabi_up_ba{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
              auto block_g_dabi_lt_ba =
                  g_dabi_lt("e,b,a,j").block(g_dabi_low_ba, g_dabi_up_ba);
              TArray G4;
              G4("e,a,i,f") = block_g_dabi_lt_ba * t2_oou_lt("f,b,j,i");
              E_O2V4_vv +=
                  TA::dot(G4("e,a,i,f"), G4("f,a,i,e") - 4.0 * G2("f,a,i,e"));
              E_O2V4_vv +=
                  TA::dot(G3("e,a,i,f"), 2.0 * G4("f,a,i,e") + G3("f,a,i,e"));
              {
                TArray G1;
                G1("e,a,i,f") = block_g_dabi_lt_ba * t2_oou_lt("f,b,i,j");
                E_O2V4_vv +=
                    TA::dot(G1("e,a,i,f"),
                            8.0 * G2("f,a,i,e") + 4.0 * G1("f,a,i,e") -
                                4.0 * G4("f,a,i,e") - 4.0 * G3("f,a,i,e"));
              }
            }
          }
        }
      }
      TA::set_default_world(global_world);
      global_world.gop.sum(E_O2V4_vv);

      {
        TArray T1;
        T1("e,i,j,f") = t2_oou_lt("e,a,j,k") * t2_oou_lt("f,a,i,k");
        TArray T3;
        T3("e,i,j,f") = t2_oou_lt("e,a,k,j") * t2_oou_lt("f,a,k,i");
        {
          TArray T4;
          T4("e,i,j,f") = t2_oou_lt("e,a,j,k") * t2_oou_lt("f,a,k,i");
          TArray G1;
          G1("e,i,j,f") = g_dabi_lt("e,a,b,i") * g_dabi_lt("f,b,a,j");
          E_O2V4_vv +=
              TA::dot(G1("e,i,j,f"), 4.0 * T1("e,i,j,f") + T3("e,i,j,f") -
                                         4.0 * T4("e,i,j,f"));
        }
        {
          TArray T2;
          T2("e,i,j,f") = t2_oou_lt("e,a,k,j") * t2_oou_lt("f,a,i,k");
          TArray G2;
          G2("e,i,j,f") = g_dabi_lt("e,a,b,i") * g_dabi_lt("f,a,b,j");
          E_O2V4_vv +=
              TA::dot(G2("e,i,j,f"), 2.0 * T2("e,i,j,f") - 2.0 * T3("e,i,j,f") -
                                         2.0 * T1("e,i,j,f"));
        }
      }

      double E_OV4 = 0;
      {
        TArray T1;
        T1("e,f") = t2_oou_lt("e,a,i,j") * t2_oou_lt("f,a,i,j");
        TArray T2;
        T2("e,f") = t2_oou_lt("e,a,i,j") * t2_oou_lt("f,a,j,i");
        {
          TArray G;
          G("e,f") = g_dabi_lt("e,a,b,i") * g_dabi_lt("f,a,b,i");
          E_OV4 = TA::dot(G("e,f"), 4.0 * T1("e,f") - 2.0 * T2("e,f"));
        }
        {
          TArray G;
          G("e,f") = g_dabi_lt("e,a,b,i") * g_dabi_lt("f,b,a,i");
          E_OV4 += TA::dot(G("e,f"), T2("e,f") - 2.0 * T1("e,f"));
        }
      }

      energy_vv = E_OV5 + E_O2V4_vv + E_OV4;

      double E_O4V = 0;
      {
        TArray T1;
        T1("m,n") = t2_ouu_lt("a,b,i,m") * t2_ouu_lt("a,b,i,n");
        TArray T2;
        T2("m,n") = t2_ouu_lt("a,b,i,m") * t2_ouu_lt("b,a,i,n");
        {
          TArray G;
          TArray T;
          G("m,n") = g_cjkl_lt("a,i,j,m") * g_cjkl_lt("a,i,j,n");
          E_O4V = TA::dot((G("m,n")), (4.0 * T1("m,n") - 2.0 * T2("m,n")));
        }
        {
          TArray G;
          TArray T;
          G("m,n") = g_cjkl_lt("a,i,j,m") * g_cjkl_lt("a,j,i,n");
          E_O4V += TA::dot((G("m,n")), (T2("m,n") - 2.0 * T1("m,n")));
        }
      }
      {
        TArray G3;
        G3("m,i,j,n") = g_cjkl_lt("a,k,j,m") * g_cjkl_lt("a,k,i,n");
        TArray G1;
        G1("m,i,j,n") = g_cjkl_lt("a,j,k,m") * g_cjkl_lt("a,i,k,n");
        {
          TArray G4;
          G4("m,i,j,n") = g_cjkl_lt("a,j,k,m") * g_cjkl_lt("a,k,i,n");
          TArray T;
          T("m,i,j,n") = t2_ouu_lt("a,b,i,m") * t2_ouu_lt("b,a,j,n");
          E_O4V += TA::dot(
              (4.0 * G1("m,i,j,n") + G3("m,i,j,n") - 4.0 * G4("m,i,j,n")),
              (T("m,i,j,n")));
        }
        {
          TArray G2;
          G2("m,i,j,n") = g_cjkl_lt("a,k,j,m") * g_cjkl_lt("a,i,k,n");
          TArray T;
          T("m,i,j,n") = t2_ouu_lt("a,b,i,m") * t2_ouu_lt("a,b,j,n");
          E_O4V += TA::dot(
              (2.0 * G2("m,i,j,n") - 2.0 * G3("m,i,j,n") - 2.0 * G1("m,i,j,n")),
              (T("m,i,j,n")));
        }
      }

      double E_O4V2 = 0.0;
      {
        TArray T1;
        T1("m,a,i,n") = t2_ouu_lt("b,a,j,m") * g_cjkl_lt("b,i,j,n");
        TArray T2;
        T2("m,a,i,n") = t2_ouu_lt("a,b,j,m") * g_cjkl_lt("b,i,j,n");
        TArray T3;
        T3("m,a,i,n") = t2_ouu_lt("b,a,j,m") * g_cjkl_lt("b,j,i,n");
        TArray T4;
        T4("m,a,i,n") = t2_ouu_lt("a,b,j,m") * g_cjkl_lt("b,j,i,n");

        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,j,i,m") * t2_ouu_lt("a,b,j,n");

          E_O4V2 =
              (TA::dot((G("m,a,i,n")), (8.0 * T1("m,a,i,n") + T4("m,a,i,n") -
                                        4.0 * T3("m,a,i,n"))));
        }
        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,i,j,m") * t2_ouu_lt("b,a,j,n");

          E_O4V2 += (TA::dot((G("m,a,i,n")), (4.0 * T1("m,a,i,n"))));
        }
        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,j,i,m") * t2_ouu_lt("b,a,j,n");

          E_O4V2 +=
              (TA::dot((G("m,a,i,n")), (2.0 * T2("m,a,i,n") + T3("m,a,i,n") -
                                        4.0 * T1("m,a,i,n"))));
        }
        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,i,j,m") * t2_ouu_lt("a,b,j,n");

          E_O4V2 +=
              (TA::dot((G("m,a,i,n")), (T2("m,a,i,n") - 4.0 * T4("m,a,i,n") -
                                        4.0 * T1("m,a,i,n"))));
        }
      }

      {
        TArray T1;
        T1("m,a,b,n") = t2_ouu_lt("c,a,i,m") * t2_ouu_lt("c,b,i,n");
        TArray T2;
        T2("m,a,b,n") = t2_ouu_lt("a,c,i,m") * t2_ouu_lt("c,b,i,n");
        TArray T3;
        T3("m,a,b,n") = t2_ouu_lt("c,a,i,m") * t2_ouu_lt("b,c,i,n");
        TArray T4;
        T4("m,a,b,n") = t2_ouu_lt("a,c,i,m") * t2_ouu_lt("b,c,i,n");

        {
          TArray G;
          G("m,a,b,n") = g_cjkl_lt("b,i,j,m") * g_cjkl_lt("a,j,i,n");

          E_O4V2 += TA::dot(
              (G("m,a,b,n")),
              (4.0 * T1("m,a,b,n") + T4("m,a,b,n") - 4.0 * T3("m,a,b,n")));
        }
        {
          TArray G;
          G("m,a,b,n") = g_cjkl_lt("b,i,j,m") * g_cjkl_lt("a,i,j,n");

          E_O4V2 += TA::dot((G("m,a,b,n")),
                            (2.0 * T2("m,a,b,n") - 2.0 * T4("m,a,b,n") -
                             2.0 * T1("m,a,b,n")));
        }
      }
      energy_oo = E_O4V + E_O4V2;

      double E_O3V3_vo = 0;
      {
        TArray T1;
        TArray T2;
        TArray T3;
        TArray T4;
        T1("e,i,a,n") = t2_oou_lt("e,b,j,i") * t2_ouu_lt("a,b,j,n");
        T2("e,i,a,n") = t2_oou_lt("e,b,i,j") * t2_ouu_lt("b,a,j,n");
        T3("e,i,a,n") = t2_oou_lt("e,b,j,i") * t2_ouu_lt("b,a,j,n");
        T4("e,i,a,n") = t2_oou_lt("e,b,i,j") * t2_ouu_lt("a,b,j,n");
        {
          TArray G1;
          TArray G2;
          G1("e,i,a,n") = g_dabi_lt("e,b,a,j") * g_cjkl_lt("b,i,j,n");
          G2("e,i,a,n") = g_dabi_lt("e,a,b,j") * g_cjkl_lt("b,j,i,n");

          E_O3V3_vo = TA::dot(3.0 * G1("e,i,a,n"), (T1("e,i,a,n")));
          E_O3V3_vo += TA::dot((G1("e,i,a,n") + G2("e,i,a,n")),
                               (T1("e,i,a,n") + 4.0 * T2("e,i,a,n") -
                                2.0 * T3("e,i,a,n") - 2.0 * T4("e,i,a,n")));
        }
        {
          TArray G3;
          TArray G4;
          G3("e,i,a,n") = g_dabi_lt("e,a,b,j") * g_cjkl_lt("b,i,j,n");
          G4("e,i,a,n") = g_dabi_lt("e,b,a,j") * g_cjkl_lt("b,j,i,n");

          E_O3V3_vo += TA::dot((G3("e,i,a,n") + G4("e,i,a,n")),
                               (T3("e,i,a,n") + T4("e,i,a,n") -
                                2.0 * T1("e,i,a,n") - 2.0 * T2("e,i,a,n")));
        }
      }

      {
        TArray T1;
        TArray T2;
        TArray T3;
        TArray T4;
        T1("e,i,j,n") = t2_oou_lt("e,a,j,k") * g_cjkl_lt("a,i,k,n");
        T2("e,i,j,n") = t2_oou_lt("e,a,k,j") * g_cjkl_lt("a,i,k,n");
        T3("e,i,j,n") = t2_oou_lt("e,a,j,k") * g_cjkl_lt("a,k,i,n");
        T4("e,i,j,n") = t2_oou_lt("e,a,k,j") * g_cjkl_lt("a,k,i,n");
        {
          TArray G;
          G("e,i,j,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("b,a,j,n");

          E_O3V3_vo += TA::dot((G("e,i,j,n")),
                               (4.0 * T1("e,i,j,n") + T4("e,i,j,n") -
                                2.0 * T3("e,i,j,n") - 2.0 * T2("e,i,j,n")));
        }
        {
          TArray G;
          G("e,i,j,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("a,b,j,n");

          E_O3V3_vo += TA::dot((G("e,i,j,n")),
                               (T2("e,i,j,n") + T3("e,i,j,n") -
                                2.0 * T4("e,i,j,n") - 2.0 * T1("e,i,j,n")));
        }
      }

      double E_O2V3_vo = 0;
      {
        TArray T1;
        TArray T2;

        T1("e,n") = t2_oou_lt("e,a,i,j") * g_cjkl_lt("a,i,j,n");
        T2("e,n") = t2_oou_lt("e,a,i,j") * g_cjkl_lt("a,j,i,n");
        {
          TArray G;
          G("e,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("a,b,i,n");
          E_O2V3_vo = TA::dot(G("e,n"), 4.0 * T1("e,n") - 2.0 * T2("e,n"));
        }
        {
          TArray G;
          G("e,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("b,a,i,n");
          E_O2V3_vo += TA::dot(G("e,n"), T2("e,n") - 2.0 * T1("e,n"));
        }
      }

      energy_vo = E_O2V4_vo + E_O2V3_vo + E_O3V3_vo;

      // singles-doubles or (E_T^(5) contribution)
      TArray g_abij_lt = g_abij_laplace_transform(
          g_abij, *this->orbital_energy(), n_occ, n_frozen, x(m));
      TArray t1_lt = t1_laplace_transform(t1, *this->orbital_energy(), n_occ,
                                          n_frozen, x(m));

      double E_O3V3_v_s = 0.0;
      {
        // This implementation is at most O3V3 scaling
        // 2G1-4G3
        {
          TArray gt1;
          gt1("f,c,j,k") = g_abij_lt("a,c,j,i") * t2_oou_lt("f,a,i,k");
          TArray ggt1;
          ggt1("k,b") = gt1("f,c,j,k") * g_dabi_lt("f,c,b,j");
          E_O3V3_v_s += TA::dot((2.0 * ggt1("k,b")), t1_lt("b,k"));

          TArray ggt2;
          ggt2("k,b") = gt1("f,c,j,k") * g_dabi_lt("f,b,c,j");
          E_O3V3_v_s += TA::dot((-4.0 * ggt2("k,b")), t1_lt("b,k"));
        }

        {
          TArray gt2;
          gt2("c,f,j,k") = g_abij_lt("c,a,j,i") * t2_oou_lt("f,a,k,i");
          TArray ggt3;
          ggt3("k,b") = gt2("c,f,j,k") * g_dabi_lt("f,b,c,j");
          E_O3V3_v_s += TA::dot((-4.0 * ggt3("k,b")), t1_lt("b,k"));
        }

        {
          TArray gt3;
          gt3("f,c,j,k") = g_abij_lt("a,c,j,i") * t2_oou_lt("f,a,k,i");
          TArray ggt4;
          ggt4("k,b") = gt3("f,c,j,k") * g_dabi_lt("f,c,b,j");
          E_O3V3_v_s += TA::dot((-4.0 * ggt4("k,b")), t1_lt("b,k"));
        }
      }

      double E_O2V3_v_s = 0.0;
      {{TArray T1;
      T1("f,i") = t1_lt("a,j") * t2_oou_lt("f,a,i,j");
      TArray G1;
      TArray G4;
      G1("f,i") = g_abij_lt("a,b,j,i") * g_dabi_lt("f,a,b,j");
      G4("f,i") = g_abij_lt("a,b,i,j") * g_dabi_lt("f,a,b,j");
      E_O2V3_v_s = TA::dot((8.0 * G1("f,i") - 4.0 * G4("f,i")), T1("f,i"));
    }
    {
      TArray T2;
      T2("f,i") = t1_lt("a,j") * t2_oou_lt("f,a,j,i");
      TArray G2;
      TArray G3;
      G2("f,i") = g_abij_lt("a,b,j,i") * g_dabi_lt("f,b,a,j");
      G3("f,i") = g_abij_lt("a,b,i,j") * g_dabi_lt("f,b,a,j");
      E_O2V3_v_s += TA::dot((2.0 * G2("f,i") - 4.0 * G3("f,i")), T2("f,i"));
    }
  }
  {
    TArray T1;
    T1("f,a") = t1_lt("b,i") * g_dabi_lt("f,b,a,i");
    TArray T2;
    T2("f,a") = t1_lt("b,i") * g_dabi_lt("f,a,b,i");

    TArray G1;
    TArray G2;
    G1("f,a") = g_abij_lt("a,b,i,j") * t2_oou_lt("f,b,i,j");
    G2("f,a") = g_abij_lt("a,b,i,j") * t2_oou_lt("f,b,j,i");

    E_O2V3_v_s += TA::dot((2.0 * G1("f,a") - G2("f,a")),
                          (4.0 * T1("f,a") - 2.0 * T2("f,a")));
  }

  {
    TArray T1;
    TArray T2;
    TArray G;
    T1("i,j,k,f") = t1_lt("a,k") * t2_oou_lt("f,a,i,j");
    T2("i,j,k,f") = t1_lt("a,k") * t2_oou_lt("f,a,j,i");
    G("i,j,k,f") = g_abij_lt("a,b,i,j") * g_dabi_lt("f,a,b,k");

    E_O3V3_v_s +=
        TA::dot((G("i,j,k,f")), (2.0 * T1("i,j,k,f") - 4.0 * T2("i,j,k,f")));
  }
  {
    TArray T1;
    TArray T2;
    TArray G3;
    T1("a,i,j,f") = t1_lt("b,j") * g_dabi_lt("f,a,b,i");
    T2("a,i,j,f") = t1_lt("b,j") * g_dabi_lt("f,b,a,i");
    G3("a,i,j,f") = g_abij_lt("b,a,k,i") * t2_oou_lt("f,b,k,j");
    {
      TArray G1;
      G1("a,i,j,f") = g_abij_lt("b,a,k,i") * t2_oou_lt("f,b,j,k");
      E_O3V3_v_s +=
          TA::dot((8.0 * G1("a,i,j,f") - 4.0 * G3("a,i,j,f")), T1("a,i,j,f"));
    }
    {
      TArray G2;
      G2("a,i,j,f") = g_abij_lt("a,b,k,i") * t2_oou_lt("f,b,j,k");
      E_O3V3_v_s +=
          TA::dot((2.0 * G2("a,i,j,f") + 2.0 * G3("a,i,j,f")), T2("a,i,j,f"));
    }
  }

  energy_v_s = E_O3V3_v_s + E_O2V3_v_s;

  double E_O3V2_o_s = 0.0;
  {
    TArray T1;
    TArray T2;
    T1("n,i") = t1_lt("a,j") * g_cjkl_lt("a,i,j,n");
    T2("n,i") = t1_lt("a,j") * g_cjkl_lt("a,j,i,n");
    {
      TArray G1;
      TArray G4;
      G1("n,i") = g_abij_lt("a,b,j,i") * t2_ouu_lt("a,b,j,n");
      G4("n,i") = g_abij_lt("a,b,i,j") * t2_ouu_lt("a,b,j,n");
      E_O3V2_o_s = TA::dot((8.0 * G1("n,i") - 4.0 * G4("n,i")), T1("n,i"));
    }
    {
      TArray G2;
      TArray G3;
      G2("n,i") = g_abij_lt("a,b,j,i") * t2_ouu_lt("b,a,j,n");
      G3("n,i") = g_abij_lt("a,b,i,j") * t2_ouu_lt("b,a,j,n");
      E_O3V2_o_s += TA::dot((2.0 * G2("n,i") - 4.0 * G3("n,i")), T2("n,i"));
    }
  }
  {
    TArray T1;
    TArray T2;
    T1("n,a") = t1_lt("b,i") * t2_ouu_lt("b,a,i,n");
    T2("n,a") = t1_lt("b,i") * t2_ouu_lt("a,b,i,n");
    {
      TArray G;
      G("n,a") = g_abij_lt("a,b,i,j") * g_cjkl_lt("b,i,j,n");
      E_O3V2_o_s += TA::dot((G("n,a")), (8.0 * T1("n,a") - 4.0 * T2("n,a")));
    }
    {
      TArray G;
      G("n,a") = g_abij_lt("a,b,i,j") * g_cjkl_lt("b,j,i,n");
      E_O3V2_o_s += TA::dot((G("n,a")), (2.0 * T2("n,a") - 4.0 * T1("n,a")));
    }
  }

  double E_O4V2_o_s = 0.0;
  {
    TArray G2;
    G2("i,j,a,n") = g_abij_lt("b,a,i,k") * g_cjkl_lt("b,k,j,n");

    {
      TArray T1;
      T1("i,j,a,n") = t1_lt("b,j") * t2_ouu_lt("a,b,i,n");
      TArray G1;
      G1("i,j,a,n") = g_abij_lt("b,a,k,i") * g_cjkl_lt("b,j,k,n");
      TArray G4;
      G4("i,j,a,n") = g_abij_lt("b,a,i,k") * g_cjkl_lt("b,j,k,n");
      E_O4V2_o_s = TA::dot(
          (8.0 * G1("i,j,a,n") + 2.0 * G2("i,j,a,n") - 4.0 * G4("i,j,a,n")),
          (T1("i,j,a,n")));
    }
    {
      TArray T2;
      T2("i,j,a,n") = t1_lt("b,j") * t2_ouu_lt("b,a,i,n");
      TArray G3;
      G3("i,j,a,n") = g_abij_lt("a,b,i,k") * g_cjkl_lt("b,j,k,n");
      E_O4V2_o_s += TA::dot((-4.0 * G3("i,j,a,n") - 4.0 * G2("i,j,a,n")),
                            (T2("i,j,a,n")));
    }
  }
  {
    TArray G;
    TArray T;
    G("i,j,k,n") = g_abij_lt("a,b,i,j") * t2_ouu_lt("a,b,k,n");
    T("i,j,k,n") = t1_lt("c,k") * g_cjkl_lt("c,j,i,n");
    E_O4V2_o_s += -4.0 * (TA::dot((G("i,j,k,n")), (T("i,j,k,n"))));
  }

  double E_O3V3_o_s = 0.0;
  {
    {
      TArray G1;
      TArray T1;
      G1("i,a,b,n") = g_abij_lt("c,a,i,j") * t2_ouu_lt("b,c,j,n");
      T1("i,a,b,n") = t1_lt("b,j") * g_cjkl_lt("a,j,i,n");
      E_O3V3_o_s += 2.0 * (TA::dot((G1("i,a,b,n")), (T1("i,a,b,n"))));
    }
    {
      TArray G2;
      TArray G3;
      TArray T2;
      G2("i,a,b,n") = g_abij_lt("a,c,i,j") * t2_ouu_lt("b,c,j,n");
      G3("i,a,b,n") = g_abij_lt("a,c,i,j") * t2_ouu_lt("c,b,j,n");
      T2("i,a,b,n") = t1_lt("b,j") * g_cjkl_lt("a,i,j,n");
      E_O3V3_o_s +=
          TA::dot((2.0 * G2("i,a,b,n") - 4.0 * G3("i,a,b,n")), (T2("i,a,b,n")));
    }
    {
      TArray gt;
      gt("i,j,k,n") = g_abij_lt("a,b,i,j") * t2_ouu_lt("a,b,k,n");
      TArray ggt;
      ggt("c,k") = gt("i,j,k,n") * g_cjkl_lt("c,i,j,n");

      E_O3V3_o_s += 2.0 * (TA::dot((ggt("c,k")), (t1_lt("c,k"))));
    }
  }

  energy_o_s = E_O3V2_o_s + E_O4V2_o_s + E_O3V3_o_s;

  triple_energy += 3.0 * w(m) *
                   (2.0 * energy_vv + 2.0 * energy_oo - 4.0 * energy_vo +
                    energy_v_s - energy_o_s);

  if (t1.world().rank() == 0) {
    std::cout << "done with quadrature: " << m + 1 << " out of "
              << n_laplace_quad_ << std::endl;
  }
}

if (t1.world().rank() == 0) {
  std::cout << "Laplace-transform-(T): " << std::endl;
  std::cout << "number of quadrature points = " << n_laplace_quad_ << std::endl;
}

triple_energy = -triple_energy / (3.0 * alpha);

return triple_energy;
}

void reblock() {
  auto &lcao_factory = this->lcao_factory();
  auto &world = lcao_factory.world();

  std::size_t b_occ = occ_block_size_;
  std::size_t b_vir = unocc_block_size_;

  std::size_t occ = this->trange1_engine()->get_occ();
  std::size_t vir = this->trange1_engine()->get_vir();
  std::size_t all = this->trange1_engine()->get_all();
  std::size_t n_frozen = this->trange1_engine()->get_nfrozen();

  TA::TiledRange1 old_occ = this->trange1_engine()->get_active_occ_tr1();
  TA::TiledRange1 old_vir = this->trange1_engine()->get_vir_tr1();

  // get occupied and virtual orbitals
  auto occ_space = lcao_factory.orbital_registry().retrieve(OrbitalIndex(L"i"));
  auto vir_space = lcao_factory.orbital_registry().retrieve(OrbitalIndex(L"a"));

  if (reblock_) {
    using TRange1Engine = ::mpqc::utility::TRange1Engine;
    auto new_tr1 =
        std::make_shared<const TRange1Engine>(occ, all, b_occ, b_vir, n_frozen);

    TA::TiledRange1 new_occ = new_tr1->get_active_occ_tr1();
    TA::TiledRange1 new_vir = new_tr1->get_vir_tr1();

    mpqc::detail::parallel_print_range_info(world, new_occ, "CCSD(T) Occ");
    mpqc::detail::parallel_print_range_info(world, new_vir, "CCSD(T) Vir");

    this->trange1_engine_ = new_tr1;

    TArray occ_convert =
        math::create_diagonal_array_from_eigen<Tile, Policy>(
            world, old_occ, new_occ, 1.0);

    TArray vir_convert =
        math::create_diagonal_array_from_eigen<Tile, Policy>(
            world, old_vir, new_vir, 1.0);

    auto new_occ_space = occ_space;
    new_occ_space("k,i") = occ_space("k,j") * occ_convert("j,i");

    auto new_vir_space = vir_space;
    new_vir_space("k,a") = vir_space("k,b") * vir_convert("b,a");

    lcao_factory.orbital_registry().clear();
    lcao_factory.orbital_registry().add(new_occ_space);
    lcao_factory.orbital_registry().add(new_vir_space);

    // get t1
    auto t1 = this->t1();
    t1("a,i") = t1("b,j") * vir_convert("b,a") * occ_convert("j,i");
    this->set_t1(t1);

    // get t2
    auto t2 = this->t2();
    t2("a,b,i,j") = t2("c,d,k,l") * vir_convert("c,a") * vir_convert("d,b") *
                    occ_convert("k,i") * occ_convert("l,j");
    this->set_t2(t2);
  }

  if (reblock_inner_) {
    auto &tr1 = this->local_trange1_engine();

    // occ inner
    tr_occ_inner_ =
        utility::compute_trange1(tr1->get_active_occ(), inner_block_size_);

    mpqc::detail::parallel_print_range_info(world, tr_occ_inner_,
                                            "CCSD(T) OCC Inner");

    auto occ_inner_convert =
        math::create_diagonal_array_from_eigen<Tile, Policy>(
            world, old_occ, tr_occ_inner_, 1.0);

    TArray inner_occ;
    inner_occ("k,i") = occ_space("k,j") * occ_inner_convert("j,i");
    OrbitalSpace<TArray> inner_occ_space =
        OrbitalSpace<TArray>(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), inner_occ);

    lcao_factory.orbital_registry().remove(OrbitalIndex(L"m"));
    lcao_factory.orbital_registry().add(inner_occ_space);

    // vir inner
    tr_vir_inner_ = utility::compute_trange1(vir, inner_block_size_);
    mpqc::detail::parallel_print_range_info(world, tr_vir_inner_,
                                            "CCSD(T) Vir Inner");
    auto vir_inner_convert =
        math::create_diagonal_array_from_eigen<Tile, Policy>(
            world, old_vir, tr_vir_inner_, 1.0);

    TArray inner_vir;
    inner_vir("k,a") = vir_space("k,b") * vir_inner_convert("b,a");
    OrbitalSpace<TArray> inner_vir_space = OrbitalSpace<TArray>(
        OrbitalIndex(L"a'"), OrbitalIndex(L"κ"), inner_vir);

    lcao_factory.orbital_registry().remove(OrbitalIndex(L"a'"));
    lcao_factory.orbital_registry().add(inner_vir_space);

    utility::print_par(world,
                       "Warning!! Using m for Inner Occupied Orbitals and a' "
                       "for Inner Virtual Orbitals! \n");
  }
}

void reblock_inner_t2(TArray &t2_left, TArray &t2_right) {
  auto &world = this->lcao_factory().world();

  auto vir_inner_convert =
      math::create_diagonal_array_from_eigen<Tile, Policy>(
          world, t2_left.trange().data()[0], tr_vir_inner_, 1.0);

  auto occ_inner_convert =
      math::create_diagonal_array_from_eigen<Tile, Policy>(
          world, t2_right.trange().data()[3], tr_occ_inner_, 1.0);

  t2_left("d,a,i,j") = t2_left("c,a,i,j") * vir_inner_convert("c,d");

  t2_right("a,b,i,l") = t2_right("a,b,i,j") * occ_inner_convert("j,l");
}

// const TArray get_Xab() {
//  TArray result;
//  TArray sqrt = this->ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
//  TArray three_center;
//  if (reblock_inner_) {
//    three_center = this->lcao_factory().compute(L"(Κ|G|a' b)");
//  } else {
//    three_center = this->lcao_factory().compute(L"(Κ|G|a b)");
//  }
//  result("K,a,b") = sqrt("K,Q") * three_center("Q,a,b");
//  return result;
//}
//
///// get three center integral (X|ai)
// const TArray get_Xai() {
//  TArray result;
//  TArray sqrt = this->ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
//  TArray three_center = this->lcao_factory().compute(L"(Κ|G|a i)");
//  result("K,a,i") = sqrt("K,Q") * three_center("Q,a,i");
//  return result;
//}
//
/// <ai|jk>
const TArray get_aijk() {
  std::wstring post_fix = this->is_df() ? L"[df]" : L"";
  TArray result;

  if (reblock_inner_) {
    result = this->lcao_factory().compute(L"<i j|G|m a>" + post_fix);
  } else {
    result = this->lcao_factory().compute(L"<i j|G|k a>" + post_fix);
  }
  result("a,i,j,k") = result("i,j,k,a");
  return result;
}

/// <ab|ci>
const TArray get_abci() {
  std::wstring post_fix = this->is_df() ? L"[df]" : L"";
  TArray result;

  if (reblock_inner_) {
    result = this->lcao_factory().compute(L"<i a'|G|b c>" + post_fix);
  } else {
    result = this->lcao_factory().compute(L"<i a|G|b c>" + post_fix);
  }
  result("a,b,c,i") = result("i,a,b,c");
  return result;
}

/// <ab|ij>
const TArray get_abij() {
  std::wstring post_fix = this->is_df() ? L"[df]" : L"";
  TArray result;
  result = this->lcao_factory().compute(L"<i j|G|a b>" + post_fix);
  result("a,b,i,j") = result("i,j,a,b");
  return result;
}

private:
struct ReduceBase {
  typedef double result_type;
  typedef Tile argument_type;

  std::shared_ptr<const Eigen::VectorXd> vec_;
  std::size_t n_occ_;
  std::size_t n_frozen_;
  std::array<std::size_t, 6> offset_;

  ReduceBase(std::shared_ptr<const Eigen::VectorXd> vec, std::size_t n_occ,
             std::size_t n_frozen, std::array<std::size_t, 6> offset)
      : vec_(std::move(vec)),
        n_occ_(n_occ),
        n_frozen_(n_frozen),
        offset_(offset) {}

  ReduceBase(ReduceBase const &) = default;

  result_type operator()() const { return 0.0; }

  result_type operator()(result_type const &t) const { return t; }

  void operator()(result_type &me, result_type const &other) const {
    me += other;
  }
};

struct CCSD_T_Reduce : public ReduceBase {
  typedef typename ReduceBase::result_type result_type;
  typedef typename ReduceBase::argument_type argument_type;

  CCSD_T_Reduce(std::shared_ptr<const Eigen::VectorXd> vec, std::size_t n_occ,
                std::size_t n_frozen, std::array<std::size_t, 6> offset)
      : ReduceBase(vec, n_occ, n_frozen, offset) {}

  CCSD_T_Reduce(CCSD_T_Reduce const &) = default;

  using ReduceBase::operator();

  void operator()(result_type &me, argument_type const &tile) const {
    auto const &ens = *this->vec_;
    std::size_t n_occ = this->n_occ_;
    auto offset_ = this->offset_;
    auto n_frozen = this->n_frozen_;

    const auto a0 = tile.range().lobound()[0];
    const auto an = tile.range().upbound()[0];
    const auto b0 = tile.range().lobound()[1];
    const auto bn = tile.range().upbound()[1];
    const auto c0 = tile.range().lobound()[2];
    const auto cn = tile.range().upbound()[2];
    const auto i0 = tile.range().lobound()[3];
    const auto in = tile.range().upbound()[3];
    const auto j0 = tile.range().lobound()[4];
    const auto jn = tile.range().upbound()[4];
    const auto k0 = tile.range().lobound()[5];
    const auto kn = tile.range().upbound()[5];

    // get the offset
    const auto a_offset = offset_[0] + n_occ;
    const auto b_offset = offset_[1] + n_occ;
    const auto c_offset = offset_[2] + n_occ;
    const auto i_offset = offset_[3] + n_frozen;
    const auto j_offset = offset_[4] + n_frozen;
    const auto k_offset = offset_[5] + n_frozen;

    auto tile_idx = 0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + a_offset];
      for (auto b = b0; b < bn; ++b) {
        const auto e_ab = e_a + ens[b + b_offset];
        for (auto c = c0; c < cn; ++c) {
          const auto e_abc = e_ab + ens[c + c_offset];
          for (auto i = i0; i < in; ++i) {
            const auto e_abci = ens[i + i_offset] - e_abc;
            for (auto j = j0; j < jn; ++j) {
              const auto e_abcij = e_abci + ens[j + j_offset];
              for (auto k = k0; k < kn; ++k, ++tile_idx) {
                const auto e_abcijk = e_abcij + ens[k + k_offset];

                me += (1.0 / e_abcijk) * tile[tile_idx];
              }
            }
          }
        }
      }
    }
  }
};  // structure CCSD_T_Reduce

struct CCSD_T_ReduceSymm : public ReduceBase {
  typedef typename ReduceBase::result_type result_type;
  typedef typename ReduceBase::argument_type argument_type;

  CCSD_T_ReduceSymm(std::shared_ptr<const Eigen::VectorXd> vec,
                    std::size_t n_occ, std::size_t n_frozen,
                    std::array<std::size_t, 6> offset)
      : ReduceBase(vec, n_occ, n_frozen, offset) {}

  CCSD_T_ReduceSymm(CCSD_T_ReduceSymm const &) = default;

  using ReduceBase::operator();

  void operator()(result_type &me, argument_type const &tile) const {
    auto const &ens = *this->vec_;
    std::size_t n_occ = this->n_occ_;
    std::size_t n_frozen = this->n_frozen_;

    // get the offset
    const auto a_offset = this->offset_[0];
    const auto b_offset = this->offset_[1];
    const auto c_offset = this->offset_[2];
    const auto i_offset = this->offset_[3];
    const auto j_offset = this->offset_[4];
    const auto k_offset = this->offset_[5];

    // compute index in the whole array
    const auto a0 = tile.range().lobound()[0] + a_offset;
    const auto an = tile.range().upbound()[0] + a_offset;
    const auto b0 = tile.range().lobound()[1] + b_offset;
    const auto bn = tile.range().upbound()[1] + b_offset;
    const auto nb = bn - b0;
    const auto c0 = tile.range().lobound()[2] + c_offset;
    const auto cn = tile.range().upbound()[2] + c_offset;
    const auto nc = cn - c0;
    const auto i0 = tile.range().lobound()[3] + i_offset;
    const auto in = tile.range().upbound()[3] + i_offset;
    const auto ni = in - i0;
    const auto j0 = tile.range().lobound()[4] + j_offset;
    const auto jn = tile.range().upbound()[4] + j_offset;
    const auto nj = jn - j0;
    const auto k0 = tile.range().lobound()[5] + k_offset;
    const auto kn = tile.range().upbound()[5] + k_offset;
    const auto nk = kn - k0;

    const auto njk = nj * nk;
    const auto nijk = ni * njk;
    const auto ncijk = nc * nijk;
    const auto nbcijk = nb * ncijk;

    typename Tile::value_type tmp = 0.0;

    // use symmetry in loop, only sum result over c <= b <= a
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      const auto aa0 = (a - a0) * nbcijk;
      for (auto b = b0; b < bn && b <= a; ++b) {
        const auto e_ab = e_a + ens[b + n_occ];
        const auto bb0 = aa0 + (b - b0) * ncijk;
        for (auto c = c0; c < cn && c <= b; ++c) {
          const auto e_abc = e_ab + ens[c + n_occ];
          const auto cc0 = bb0 + (c - c0) * nijk;
          bool none_equal = (a != b && a != c && b != c);
          bool diagonal = (a == b && b == c);
          for (auto i = i0; i < in; ++i) {
            const auto e_abci = ens[i + n_frozen] - e_abc;
            const auto ii0 = cc0 + (i - i0) * njk;
            for (auto j = j0; j < jn; ++j) {
              const auto e_abcij = e_abci + ens[j + n_frozen];
              const auto jj0 = ii0 + (j - j0) * nk;
              for (auto k = k0; k < kn; ++k) {
                const auto e_abcijk = e_abcij + ens[k + n_frozen];
                const auto tile_idx = jj0 + (k - k0);

                tmp = (1.0 / e_abcijk) * tile[tile_idx];
                // 6 fold symmetry if none in a,b,c equal
                if (none_equal) {
                  tmp = 2.0 * tmp;
                }
                // if diagonal, a==b==c no symmetry
                else if (diagonal) {
                  tmp = 0;
                }
                // three fold symmetry if two in a,b,c equal
                else {
                  tmp = tmp;
                }
                me += tmp;
              }
            }
          }
        }
      }
    }
  }
};  // structure CCSD_T_ReduceSymm

};  // end of class CCSD_T

#if TA_DEFAULT_POLICY == 0
extern template class CCSD_T<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class CCSD_T<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_T_H_
