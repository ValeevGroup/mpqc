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
  bool replicate_;

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

 public:
  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords : all keywords for CCSD
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | reblock_occ | int | none | block size to reblock occ |
   * | reblock_unocc | int | none | block size to reblock unocc |
   * | reblock_inner | int | size of number of orbital | block size to reblock inner dimension, set to 0 to disable reblock inner |
   * | approach | string | coarse | coarse grain, fine grain or straight approach |
   * | replicate | bool | false | valid only with approach=coarse, if replicate integral g_cijk(smallest integral in (T)) |
   * | increase | int | 2 | valid only with approach=fine, number of block in virtual dimension to load at each virtual loop |
   * | quadrature_points | int | 4(need to benchmark it) | number of quadrature points for the Laplace transform |
   */
  // clang-format on

  CCSD_T(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    reblock_ = kv.exists("reblock_occ") || kv.exists("reblock_unocc");

    occ_block_size_ = kv.value<int>("reblock_occ", 8);
    unocc_block_size_ = kv.value<int>("reblock_unocc", 8);
    replicate_ = kv.value<bool>("replicate", false);

    // default value is size of total number of orbitals, which makes it 1 block
    std::size_t n_obs = this->wfn_world()
                            ->basis_registry()
                            ->retrieve(OrbitalIndex(L"κ"))
                            ->nfunctions();
    inner_block_size_ = kv.value<int>("reblock_inner", n_obs);
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
    this->lcao_factory().registry().purge(world);

    if(approach_ != "laplace"){
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
      ExEnv::out0() << "(T) Time in CCSD(T): " << duration0 << std::endl;
    }
  }

 private:
  double compute_ccsd_t_coarse_grain(TArray &t1, TArray &t2) {
    auto &global_world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    // get integral
    TArray g_cjkl_global = get_aijk();
    TArray g_abij = get_abij();
    TArray g_dabi = get_abci();

    // T2
    TArray t2_left = t2;
    TArray t2_right = t2;

    // if reblock_inner, needs two copy of T2 with different blocking
    if (reblock_inner_) {
      reblock_inner_t2(t2_left, t2_right);
    }

    // get trange1
    auto tr_occ = this->trange1_engine_->get_active_occ_tr1();
    auto tr_vir = this->trange1_engine_->get_vir_tr1();

    auto n_tr_occ = this->trange1_engine_->get_active_occ_blocks();
    auto n_tr_vir = this->trange1_engine_->get_vir_blocks();

    // TiledRange1 for occ, unocc and inner contraction space
    auto n_tr_occ_inner = n_tr_occ;
    auto n_tr_vir_inner = n_tr_vir;
    std::cout << n_tr_vir_inner << std::endl;
    if (reblock_inner_) {
      n_tr_occ_inner = tr_occ_inner_.tiles_range().second;
      n_tr_vir_inner = tr_vir_inner_.tiles_range().second;
      std::cout << n_tr_vir_inner << std::endl;
    }

    std::size_t vir_block_size = this->trange1_engine_->get_vir_block_size();
    std::size_t n_occ = this->trange1_engine_->get_active_occ();
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
      world_ptr = std::shared_ptr<madness::World>(new madness::World(comm));
      tmp_ptr = world_ptr.get();
    } else {
      tmp_ptr = &global_world;
    }
    auto &this_world = *tmp_ptr;
    global_world.gop.fence();

    TArray g_cjkl(
        this_world, g_cjkl_global.trange(), g_cjkl_global.shape(),
        Policy::default_pmap(this_world,
                             g_cjkl_global.trange().tiles_range().volume()));
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
    if (replicate_ && size > 1) {
      // replicate g_cjkl
      g_cjkl("c,j,k,l") = g_cjkl_global("c,j,k,l").set_world(this_world);
    } else {
      // don't replicate g_cjkl
      g_cjkl = g_cjkl_global;
    }

    // lambda function to compute t3
    auto compute_t3 = [&](std::size_t a, std::size_t b, std::size_t c,
                          TArray &t3) {
      // index
      std::size_t a_low = a;
      std::size_t a_up = a + 1;
      std::size_t b_low = b;
      std::size_t b_up = b + 1;
      std::size_t c_low = c;
      std::size_t c_up = c + 1;
      typedef std::vector<std::size_t> block;

      {
        // block for g_dbai
        block g_dabi_low{0, a_low, b_low, 0};
        block g_dabi_up{n_tr_vir_inner, a_up, b_up, n_tr_occ};

        std::cout << "a_low = " << a_low << std::endl;
        std::cout << "g_dabi.trange() = " << g_dabi.trange() << std::endl;
        auto block_g_dabi = g_dabi("d,a,b,i").block(g_dabi_low, g_dabi_up);

        // block for t2_cdk
        block t2_dcjk_low{0, c_low, 0, 0};
        block t2_dcjk_up{n_tr_vir_inner, c_up, n_tr_occ, n_tr_occ};

        auto block_t2_dcjk = t2_left("d,c,j,k").block(t2_dcjk_low, t2_dcjk_up);

        std::cout << "t2_left.trange() = " << t2_left.trange() << std::endl;

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
              ((block_g_dabi * block_t2_dcjk).set_world(this_world) -
               (block_t2_abil * block_g_cjkl).set_world(this_world))
                  .set_world(this_world);
        } else {
          t3("a,b,i,c,j,k") =
              ((block_g_dabi * block_t2_dcjk).set_world(this_world) -
               (block_t2_abil * block_g_cjkl).set_world(this_world))
                  .set_world(this_world);
        }
        std::cout << "t3.trange() = " << t3.trange() << std::endl;
      }
    };

    // lambda function to compute v3
    auto compute_v3 = [&](std::size_t a, std::size_t b, std::size_t c,
                          TArray &v3) {
      // index
      std::size_t a_low = a;
      std::size_t a_up = a + 1;
      std::size_t b_low = b;
      std::size_t b_up = b + 1;
      std::size_t c_low = c;
      std::size_t c_up = c + 1;
      typedef std::vector<std::size_t> block;

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
          v3("b,c,j,k,a,i") +=
              (block_g_bcjk * block_t_ai).set_world(this_world);
        } else {
          v3("b,c,j,k,a,i") = (block_g_bcjk * block_t_ai).set_world(this_world);
        }
      }
    };

    // start (T) calculation
    double triple_energy = 0.0;

    // time spend in reduction
    double reduce_time = 0.0;
    // time spend in contraction
    double contraction_time = 0.0;
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

    // start loop over a, b, c
    for (auto a = 0; a < n_tr_vir; ++a) {
      for (auto b = 0; b <= a; ++b) {
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
          compute_t3(a, b, c, t3);
          time01 = mpqc::now(this_world, accurate_time);
          contraction_time += mpqc::duration_in_s(time00, time01);

          // acbikj contribution
          // g^{da}_{ci}*t^{db}_{kj} - g^{bk}_{jl}*t^{ac}_{il}
          {
            t3("a,c,i,b,k,j") = t3("a,b,i,c,j,k");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_t3(a, c, b, t3);
            time01 = mpqc::now(this_world, accurate_time);
            contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // cabkij contribution
          // g^{dc}_{ak}*t^{db}_{ij} - g^{bi}_{jl}*t^{ca}_{kl}
          {
            t3("c,a,k,b,i,j") = t3("a,c,i,b,k,j");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_t3(c, a, b, t3);
            time01 = mpqc::now(this_world, accurate_time);
            contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // cbakji contribution
          // g^{dc}_{bk}*t^{da}_{ji} - g^{aj}_{il}*t^{cb}_{kl}
          {
            t3("c,b,k,a,j,i") = t3("c,a,k,b,i,j");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_t3(c, b, a, t3);
            time01 = mpqc::now(this_world, accurate_time);
            contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // bcajki contribution
          // g^{db}_{cj}*t^{da}_{ki} - g^{ak}_{il}*t^{bc}_{jl}
          {
            t3("b,c,j,a,k,i") = t3("c,b,k,a,j,i");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_t3(b, c, a, t3);
            time01 = mpqc::now(this_world, accurate_time);
            contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // bacjik contribution
          // g^{db}_{aj}*t^{dc}_{ik} - g^{ci}_{kl}*t^{ba}_{jl}
          {
            t3("b,a,j,c,i,k") = t3("b,c,j,a,k,i");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_t3(b, a, c, t3);
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
            contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // acbikj contribution
          // g^{ac}_{ik}*t^{b}_{j}
          {
            v3("a,c,i,k,b,j") = v3("b,c,j,k,a,i");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_v3(b, a, c, v3);
            time01 = mpqc::now(this_world, accurate_time);
            contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // abcijk contribution
          // g^{ab}_{ij}*t^{c}_{k}
          {
            v3("a,b,i,j,c,k") = v3("a,c,i,k,b,j");
            time00 = mpqc::now(this_world, accurate_time);
            permutation_time += mpqc::duration_in_s(time01, time00);

            compute_v3(c, a, b, v3);
            time01 = mpqc::now(this_world, accurate_time);
            contraction_time += mpqc::duration_in_s(time00, time01);
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
                this->orbital_energy_, this->trange1_engine_->get_occ(),
                this->trange1_engine_->get_nfrozen(), offset);
            tmp_energy = result("a,b,c,i,j,k").reduce(ccsd_t_reduce);
            time01 = mpqc::now(this_world, accurate_time);
            reduce_time += mpqc::duration_in_s(time00, time01);
            tmp_energy *= 2;
          }
          // boundary condition
          else {
            time00 = mpqc::now(this_world, accurate_time);
            auto ccsd_t_reduce = CCSD_T_ReduceSymm(
                this->orbital_energy_, this->trange1_engine_->get_occ(),
                this->trange1_engine_->get_nfrozen(), offset);
            tmp_energy = result("a,b,c,i,j,k").reduce(ccsd_t_reduce);
            time01 = mpqc::now(this_world, accurate_time);
            reduce_time += mpqc::duration_in_s(time00, time01);
          }

          triple_energy += tmp_energy;
        }  // loop of c
      }    // loop of b

      // print the progress
      if (rank == 0) {
        print_progress(a, a + 1, progress_points);
      }
    }  // loop of a
    this_world.gop.fence();
    global_world.gop.fence();

    if (this->print_detail()) {
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
          std::cout << "Reduce Time: " << reduce_time << " S" << std::endl
                    << std::endl;
        }
      }
    }

    // print out all process time
    global_world.gop.sum(iter);
    global_world.gop.sum(permutation_time);
    global_world.gop.sum(contraction_time);
    global_world.gop.sum(reduce_time);

    ExEnv::out0() << "Process All Time: " << std::endl;
    ExEnv::out0() << "Iter: " << iter << std::endl;
    ExEnv::out0() << "Permutation Time: " << permutation_time << " S"
                  << std::endl;
    ExEnv::out0() << "Contraction Time: " << contraction_time << " S"
                  << std::endl;
    ExEnv::out0() << "Reduce Time: " << reduce_time << " S" << std::endl
                  << std::endl;

    global_world.gop.sum(triple_energy);

    //    ExEnv::out0() << "(T) Energy: " << triple_energy << std::endl
    //                  << std::endl;

    // manually clean replicated array
    if (size > 1) {
      t1_this = TArray();
      g_cjkl = TArray();

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

    // get trange1
    auto tr_occ = this->trange1_engine_->get_active_occ_tr1();
    auto tr_vir = this->trange1_engine_->get_vir_tr1();

    auto n_tr_occ = this->trange1_engine_->get_active_occ_blocks();
    auto n_tr_vir = this->trange1_engine_->get_vir_blocks();
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

    std::size_t occ_block_size = this->trange1_engine_->get_occ_block_size();
    std::size_t vir_block_size = this->trange1_engine_->get_vir_block_size();
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
                this->orbital_energy_, this->trange1_engine_->get_occ(),
                this->trange1_engine_->get_nfrozen(), offset);
            tmp_energy = ((t3("a,b,c,i,j,k") + v3("a,b,c,i,j,k")) *
                          (4.0 * t3("a,b,c,i,j,k") + t3("a,b,c,k,i,j") +
                           t3("a,b,c,j,k,i") -
                           2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") +
                                t3("a,b,c,j,i,k"))))
                             .reduce(ccsd_t_reduce);

            tmp_energy *= 2;
          } else {
            auto ccsd_t_reduce = CCSD_T_ReduceSymm(
                this->orbital_energy_, this->trange1_engine_->get_occ(),
                this->trange1_engine_->get_nfrozen(), offset);
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
        print_progress(a, a + increase, progress_points);
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

  // compute and store all t3 amplitudes, not recommanded for performance
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

    auto ccsd_t_reduce =
        CCSD_T_Reduce(this->orbital_energy_, this->trange1_engine_->get_occ(),
                      this->trange1_engine_->get_nfrozen(), offset);

    double triple_energy =
        ((t3("a,b,c,i,j,k") + v3("a,b,c,i,j,k")) *
         (4.0 * t3("a,b,c,i,j,k") + t3("a,b,c,k,i,j") + t3("a,b,c,j,k,i") -
          2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") + t3("a,b,c,j,i,k"))))
            .reduce(ccsd_t_reduce);
    triple_energy = triple_energy / 3.0;
    return triple_energy;
  }

  // this function uses orthogonal polynomial approach for obtaining quadrature
  // roots and weights.
  // Taken from paper: Gaussian Quadrature and the Eigenvalue Problem by John A.
  // Gubner.
  // http://gubner.ece.wisc.edu/gaussquad.pdf
  // The code is outlined at Example 15: Legendre polynomials
  // void gauss_quad(int N, double a, double b, Eigen::VectorXd &w,
  //                 Eigen::VectorXd &x) {
  //   Eigen::MatrixXd J;
  //   J.setZero(N, N);
  //   for (auto i = 0; i < N; i++) {
  //     if (i < N - 1) {
  //       J(i, i + 1) = sqrt(1 / (4 - pow(i + 1, -2)));
  //     }
  //   }
  //   Eigen::MatrixXd Jfin = J + J.transpose();

  //   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Jfin);
  //   x = es.eigenvalues();
  //   Eigen::MatrixXd V = es.eigenvectors();

  //   for (auto i = 0; i < N; i++) {
  //     w(i) = 0.5 * 2.0 * (b - a) * V(0, i) * V(0, i);
  //     x(i) = (b - a) * 0.5 * x(i) + (b + a) * 0.5;
  //   }
  // }

  // performs Laplace transform perturbative triple correction to CCSD energy
  double compute_ccsd_t_laplace_transform(const TArray &t1, const TArray &t2) {

    auto &world = this->wfn_world()->world();
    bool accurate_time = this->lcao_factory().accurate_time();

    // get integral
    TArray g_cjkl = get_aijk();
    TArray g_dabi = get_abci();
    TArray g_abij = get_abij();


    //obtaining DF-integrals
    TArray Xab;
    TArray Xai;
    if (this->is_df()) {
      std::cout << "DF calculation" << std::endl;
      /// get three center integral (X|ab)
      //TArray Xab;
      TArray sqrt =
          this->lcao_factory().ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
      TArray Kab = this->lcao_factory().compute(L"(Κ|G|a b)");
      Xab("K,a,b") = sqrt("K,Q") * Kab("Q,a,b");

      //TArray Xai;
      TArray Kai = this->lcao_factory().compute(L"(Κ|G|a i)");
      Xai("K,a,i") = sqrt("K,Q") * Kai("Q,a,i");
    }


    // definition of orbital spaces
    auto n_occ = this->trange1_engine()->get_occ();
    auto n_frozen = this->trange1_engine()->get_nfrozen();

    // copying orbital energies into Eigen vector
    Eigen::VectorXd &e_orb = *this->orbital_energy();

    // definition of alpha parameter required for the Gaussian quadrature
    // defined as in paper: Constans, Pere, Philippe Y. Ayala, and Gustavo E.
    // Scuseria.
    //"Scaling reduction of the perturbative triples correction (T) to coupled
    //cluster
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
    // 0 and 1 are integration intervals
    // gauss_quad(n, 0, 1, w, x);
    mpqc::math::gauss_legendre(n, w, x);

    double triple_energy = 0.0;
    double triple_energys = 0.0;
    double triple_energyo = 0.0;
    double triple_energyos = 0.0;
    double triple_energyuo = 0.0;
    double triple_energy_total = 0.0;

    mpqc::time_point time00;
    mpqc::time_point time01;
    mpqc::time_point time02;
    mpqc::time_point time03;
    mpqc::time_point time04;
    mpqc::time_point time05;
    mpqc::time_point time06;
    mpqc::time_point time07;
    mpqc::time_point time08;
    mpqc::time_point time09;
    mpqc::time_point time10;
    mpqc::time_point time11;
    mpqc::time_point time12;
    mpqc::time_point time13;
    mpqc::time_point time14;
    mpqc::time_point time15;
    mpqc::time_point time16;
    mpqc::time_point time17;
    mpqc::time_point time18;
    mpqc::time_point time19;
    mpqc::time_point time20;
    mpqc::time_point time21;
    mpqc::time_point time22;
    mpqc::time_point time23;
    mpqc::time_point time24;
    mpqc::time_point time25;
    mpqc::time_point time26;
    mpqc::time_point time27;
    mpqc::time_point time28;
    mpqc::time_point time29;
    mpqc::time_point time30;
    mpqc::time_point time31;
    mpqc::time_point time32;
    mpqc::time_point time33;
    mpqc::time_point time34;
    mpqc::time_point time35;
    mpqc::time_point time36;
    mpqc::time_point time37;

    double int_transform1 = 0.0;
    double int_transform2 = 0.0;
    double int_transform3 = 0.0;
    double time_t2_vv = 0.0;
    double time_t2_ov = 0.0;
    double time_t2_oo = 0.0;
    double time_t1 = 0.0;

    double time_T_OV5 = 0.0;
    double time_G_OV5 = 0.0;
    double time_trace = 0.0;
    double time_gg1 = 0.0;
    double time_gg2 = 0.0;
    double time_ggg1 = 0.0;
    double time_ggg2 = 0.0;
    double time_G1 = 0.0;
    double time_G2 = 0.0;


    // get trange1
    auto n_tr_occ = this->trange1_engine_->get_active_occ_blocks();
    auto n_tr_vir = this->trange1_engine_->get_vir_blocks();


    auto& global_world = this->wfn_world()->world();
    // split global_world
    const auto rank = global_world.rank();
    const auto size = global_world.size();

    madness::World *tmp_ptr;
    std::shared_ptr<madness::World> world_ptr;

    if (size > 1) {
      SafeMPI::Group group = global_world.mpi.comm().Get_group().Incl(1, &rank);
      SafeMPI::Intracomm comm = global_world.mpi.comm().Create(group);
      world_ptr = std::shared_ptr<madness::World>(new madness::World(comm));
      tmp_ptr = world_ptr.get();
    } else {
      tmp_ptr = &global_world;
    }
    auto &this_world = *tmp_ptr;
    global_world.gop.fence();

    // loop over number of quadrature points
    for (auto m = 0; m < n; m++) {

      time00 = mpqc::now(world, accurate_time);
      TArray g_dabi_lt = g_dabi_laplace_transform(
          g_dabi, *this->orbital_energy(), n_occ, n_frozen, x(m));
      TArray t2_oou_lt = t2_oou_laplace_transform(t2, *this->orbital_energy(),
                                                  n_occ, n_frozen, x(m));
      time01 = mpqc::now(world, accurate_time);
      int_transform1 += mpqc::duration_in_s(time00, time01);


      time04 = mpqc::now(world, accurate_time);
      TArray g_cjkl_lt = g_cjkl_laplace_transform(
          g_cjkl, *this->orbital_energy(), n_occ, n_frozen, x(m));
      TArray t2_ouu_lt = t2_ouu_laplace_transform(t2, *this->orbital_energy(),
          n_occ, n_frozen, x(m));
      time05 = mpqc::now(world, accurate_time);
      int_transform2 += mpqc::duration_in_s(time04, time05);


      time10 = mpqc::now(world, accurate_time);
      TArray g_abij_lt = g_abij_laplace_transform(
          g_abij, *this->orbital_energy(), n_occ, n_frozen, x(m));
      TArray t1_lt = t1_laplace_transform(t1, *this->orbital_energy(), n_occ,
          n_frozen, x(m));
      time11 = mpqc::now(world, accurate_time);
      int_transform3 += mpqc::duration_in_s(time10, time11);

      TArray Xab_lt;
      TArray Xai_lt;
      if (this->is_df()) {
        Xab_lt = Xab_laplace_transform(
            Xab, *this->orbital_energy(), n_occ, n_frozen, x(m));

        Xai_lt = Xai_laplace_transform(
            Xai, *this->orbital_energy(), n_occ, n_frozen, x(m));

      }

      this->wfn_world()->world().gop.fence();
      double energy_m = 0.0;
      double energy_ms = 0.0;
      double energy_mo = 0.0;
      double energy_mos = 0.0;
      double energy_muo = 0.0;
      double Wijkabc = 0.0;

      time02 = mpqc::now(world, accurate_time);

      //computation of the OV5 terms
      double E_O2V4_vo = 0;
      /*double E_OV5 = 0;
      {
        TArray T1;
        TArray T2;

        T1("e,b,a,f") = t2_oou_lt("e,a,i,j")*t2_oou_lt("f,b,i,j");
        T2("e,b,a,f") = t2_oou_lt("e,a,i,j")*t2_oou_lt("f,b,j,i");
        {
          TArray G;
          TArray g_tmp;
          g_tmp("e,c,b,i") = g_dabi_lt("e,c,b,i") - g_dabi_lt("e,b,c,i");
          G("e,b,a,f") = g_tmp("e,c,b,i")*g_dabi_lt("f,c,a,i");

          E_OV5 = TA::dot(G("e,b,a,f"),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));
        }
        //with DF
        {
          TArray gg;
          gg("X,Y") = Xai_lt("X,c,i") * Xai_lt("Y,c,i");

          TArray ggg;
          ggg("e,b,Y") = gg("X,Y") * Xab_lt("X,e,b");

          TArray G1;
          G1("e,b,a,f") = ggg("e,b,Y") * Xab_lt("Y,f,a");

          E_OV5 = TA::dot(G1("e,b,a,f"),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));

          TArray gg1;
          gg1("e,Y,i,X") = Xab_lt("X,e,c") * Xai_lt("Y,c,i");

          TArray ggg1;
          ggg1("e,b,Y") = gg1("e,Y,i,X") * Xai_lt("X,b,i");

          TArray G2;
          G2("e,b,a,f") = ggg1("e,b,Y") * Xab_lt("Y,f,a");

          E_OV5 -= TA::dot(G2("e,b,a,f"),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));
        }

        {
          TArray G;
          G("e,b,a,f") = g_dabi_lt("e,b,c,i")*g_dabi_lt("f,a,c,i");
          E_OV5 += TA::dot(G("e,b,a,f"),(T2("e,b,a,f") - 2.0*T1("e,b,a,f")));

        }
        std::cout << "E_OV5 = " << E_OV5 << std::endl;
      }*/

      //double E_OV5_block = 0.0;
      double E_OV5 = 0;
      TA::set_default_world(this_world);
      std::size_t global_iter = 0;
      for (auto e = 0; e < n_tr_vir; ++e) {
        for (auto f = 0; f < n_tr_vir; ++f) {

          ++ global_iter;

          if( global_iter % size != rank) continue;

          std::size_t e_low = e;
          std::size_t e_up = e + 1;
          std::size_t f_low = f;
          std::size_t f_up = f + 1;

          typedef std::vector<std::size_t> block;

          //blocking g_dabi over e
          block g_dabi_low_e{e_low, 0, 0, 0};
          block g_dabi_up_e{e_up, n_tr_vir, n_tr_vir, n_tr_occ};

          auto block_g_dabi_lt_e = g_dabi_lt("e,b,c,i").block(g_dabi_low_e, g_dabi_up_e);

          //blocking g_dabi over f
          block g_dabi_low_f{f_low, 0, 0, 0};
          block g_dabi_up_f{f_up, n_tr_vir, n_tr_vir, n_tr_occ};


          //blocking t2 over e (required for both T1 and T2 intermediates)
          block t2_oou_lt_low_e{e_low, 0, 0, 0};
          block t2_oou_lt_up_e{e_up, n_tr_vir, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_e_T12 = t2_oou_lt("e,a,i,j").block(t2_oou_lt_low_e, t2_oou_lt_up_e);

          //blocking t2 over e (required for T2 intermediate)
          block t2_oou_lt_low_f{f_low, 0, 0, 0};
          block t2_oou_lt_up_f{f_up, n_tr_vir, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_f_T2 = t2_oou_lt("f,b,j,i").block(t2_oou_lt_low_f, t2_oou_lt_up_f);

          //blocking t2 over e (required for T1 intermediate)
          auto block_t2_oou_lt_f_T1 = t2_oou_lt("f,b,i,j").block(t2_oou_lt_low_f, t2_oou_lt_up_f);

          TArray T1;
          TArray T2;
          time20 = mpqc::now(world, accurate_time);
          T1("e,b,a,f") = block_t2_oou_lt_e_T12*block_t2_oou_lt_f_T1;
          T2("e,b,a,f") = block_t2_oou_lt_e_T12*block_t2_oou_lt_f_T2;
          time21 = mpqc::now(world, accurate_time);
          time_T_OV5 += mpqc::duration_in_s(time20, time21);
          {
            TArray G;

            time22 = mpqc::now(world, accurate_time);
            auto block_g_dabi_lt_f = g_dabi_lt("f,a,c,i").block(g_dabi_low_f, g_dabi_up_f);


            G("e,b,a,f") = block_g_dabi_lt_e*block_g_dabi_lt_f;
            time23 = mpqc::now(world, accurate_time);
            time_G_OV5 += mpqc::duration_in_s(time22, time23);

            E_OV5 += TA::dot((G("e,b,a,f")),(T2("e,b,a,f") - 2.0*T1("e,b,a,f")));
          }
          {
            //if (this->is_df()) {

              time22 = mpqc::now(world, accurate_time);
              auto n_tr_aux = Xab.range().upbound()[0];

              block Xab_low_e{0, e_low, 0};
              block Xab_up_e{n_tr_aux, e_up, n_tr_vir};
              auto block_Xab_lt_e = Xab_lt("X,e,b").block(Xab_low_e, Xab_up_e);

              block Xab_low_f{0, f_low, 0};
              block Xab_up_f{n_tr_aux, f_up, n_tr_vir};
              auto block_Xab_lt_f = Xab_lt("Y,f,a").block(Xab_low_f, Xab_up_f);

              time25 = mpqc::now(world, accurate_time);
              TArray gg1;
              gg1("X,Y") = Xai_lt("X,c,i") * Xai_lt("Y,c,i");
              time26 = mpqc::now(world, accurate_time);
              time_gg1 += mpqc::duration_in_s(time25, time26);

              time27 = mpqc::now(world, accurate_time);
              TArray ggg1;
              ggg1("e,Y,b") = gg1("X,Y") * block_Xab_lt_e;
              time28 = mpqc::now(world, accurate_time);
              time_ggg1 += mpqc::duration_in_s(time27, time28);

              auto block_Xab_lt_ec = Xab_lt("X,e,c").block(Xab_low_e, Xab_up_e);

              time29 = mpqc::now(world, accurate_time);
              TArray gg2;
              gg2("e,Y,i,X") = block_Xab_lt_ec * Xai_lt("Y,c,i");
              time30 = mpqc::now(world, accurate_time);
              time_gg2 += mpqc::duration_in_s(time29, time30);

              time31 = mpqc::now(world, accurate_time);
              TArray ggg2;
              ggg2("e,Y,b") = gg2("e,Y,i,X") * Xai_lt("X,b,i");
              time32 = mpqc::now(world, accurate_time);
              time_ggg2 += mpqc::duration_in_s(time31, time32);

              //TArray G;
              //G("e,b,a,f") = ggg1("e,Y,b") * block_Xab_lt_f - ggg2("e,Y,b") * block_Xab_lt_f;
              time33 = mpqc::now(world, accurate_time);
              TArray G1;
              G1("e,b,a,f") = ggg1("e,Y,b") * block_Xab_lt_f;
              time34 = mpqc::now(world, accurate_time);
              time_G1 += mpqc::duration_in_s(time33, time34);

              time35 = mpqc::now(world, accurate_time);
              TArray G2;
              G2("e,b,a,f") = ggg2("e,Y,b") * block_Xab_lt_f;
              time36 = mpqc::now(world, accurate_time);
              time_G2 += mpqc::duration_in_s(time35, time36);
              //E_OV5 += TA::dot((G("e,b,a,f")),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));

              E_OV5 += TA::dot((G1("e,b,a,f") - G2("e,b,a,f")),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));

            /*} else {
              TArray G;

              time22 = mpqc::now(world, accurate_time);
              auto block_g_dabi_lt_f = g_dabi_lt("f,c,a,i").block(g_dabi_low_f, g_dabi_up_f);
              auto block_g_dabi_lt_e1 = g_dabi_lt("e,c,b,i").block(g_dabi_low_e, g_dabi_up_e);

              G("e,b,a,f") = block_g_dabi_lt_e1*block_g_dabi_lt_f - block_g_dabi_lt_e*block_g_dabi_lt_f;
              time23 = mpqc::now(world, accurate_time);
              time_G_OV5 += mpqc::duration_in_s(time22, time23);
              E_OV5 += TA::dot((G("e,b,a,f")),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));
            }*/
          }
        }
      }
      TA::set_default_world(global_world);
      global_world.gop.sum(E_OV5);

      //double E_OV5 = 0;
      double E_O2V4_2 = 0;
      double E_O2V4_S = 0.0;
      double E_O3V3_S2 = 0.0;
      TA::set_default_world(this_world);
      global_iter = 0;
      for (auto a = 0; a < n_tr_vir; ++a) {
        for (auto b = 0; b < n_tr_vir; ++b) {
          ++ global_iter;

          if( global_iter % size != rank) continue;

          std::size_t a_low = a;
          std::size_t a_up = a + 1;
          std::size_t b_low = b;
          std::size_t b_up = b + 1;

          typedef std::vector<std::size_t> block;

          //blocking g_dabi over b
          block g_dabi_low_b{0, b_low, 0, 0};
          block g_dabi_up_b{n_tr_vir, b_up, n_tr_vir, n_tr_occ};

          auto block_g_dabi_lt_b = g_dabi_lt("e,b,c,i").block(g_dabi_low_b, g_dabi_up_b);

          //blocking g_dabi over a
          block g_dabi_low_a{0, a_low, 0, 0};
          block g_dabi_up_a{n_tr_vir, a_up, n_tr_vir, n_tr_occ};
          auto block_g_dabi_lt_a = g_dabi_lt("f,a,c,i").block(g_dabi_low_a, g_dabi_up_a);

          //blocking t2 over a (required for both T1 and T2 intermediates)
          block t2_oou_lt_low_a{0, a_low, 0, 0};
          block t2_oou_lt_up_a{n_tr_vir, a_up, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_a_T12 = t2_oou_lt("e,a,i,j").block(t2_oou_lt_low_a, t2_oou_lt_up_a);

          //blocking t2 over b (required for T2 intermediate)
          block t2_oou_lt_low_b{0, b_low, 0, 0};
          block t2_oou_lt_up_b{n_tr_vir, b_up, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_b_T2 = t2_oou_lt("f,b,j,i").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

          //blocking t2 over b (required for T1 intermediate)
          auto block_t2_oou_lt_b_T1 = t2_oou_lt("f,b,i,j").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

          /*TArray T1;
          TArray T2;

          time20 = mpqc::now(world, accurate_time);


          T1("e,b,a,f") = block_t2_oou_lt_a_T12*block_t2_oou_lt_b_T1;
          T2("e,b,a,f") = block_t2_oou_lt_a_T12*block_t2_oou_lt_b_T2;

          time21 = mpqc::now(world, accurate_time);
          time_T_OV5 += mpqc::duration_in_s(time20, time21);

          {
            TArray G;

            time22 = mpqc::now(world, accurate_time);
            G("e,b,a,f") = block_g_dabi_lt_b*block_g_dabi_lt_a;
            time23 = mpqc::now(world, accurate_time);
            time_G_OV5 += mpqc::duration_in_s(time22, time23);

            time24 = mpqc::now(world, accurate_time);
            E_OV5 += TA::dot((G("e,b,a,f")),(T2("e,b,a,f") - 2.0*T1("e,b,a,f")));
            time25 = mpqc::now(world, accurate_time);
            time_trace += mpqc::duration_in_s(time24, time25);
          }
          {
            TArray G;

            time22 = mpqc::now(world, accurate_time);
            block g_dabi_low_ca{0, 0, a_low, 0};
            block g_dabi_up_ca{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
            auto block_g_dabi_lt_ca = g_dabi_lt("f,c,a,i").block(g_dabi_low_ca, g_dabi_up_ca);

            block g_dabi_low_cb{0, 0, b_low, 0};
            block g_dabi_up_cb{n_tr_vir, n_tr_vir, b_up, n_tr_occ};
            auto block_g_dabi_lt_cb = g_dabi_lt("e,c,b,i").block(g_dabi_low_cb, g_dabi_up_cb);

            G("e,b,a,f") = block_g_dabi_lt_cb*block_g_dabi_lt_ca - block_g_dabi_lt_b*block_g_dabi_lt_ca;
            time23 = mpqc::now(world, accurate_time);
            time_G_OV5 += mpqc::duration_in_s(time22, time23);

            time24 = mpqc::now(world, accurate_time);
            E_OV5 += TA::dot((G("e,b,a,f")),(4.0*T2("e,b,a,f") - 2.0*T1("e,b,a,f")));
            time25 = mpqc::now(world, accurate_time);
            time_trace += mpqc::duration_in_s(time24, time25);
          }*/

          //Mixed term contributions
          {
            TArray T1;
            TArray T2;

            auto block_t2_oou_lt_eb = t2_oou_lt("e,b,i,j").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

            block g_cjkl_low_a{a_low, 0, 0, 0};
            block g_cjkl_up_a{a_up, n_tr_occ, n_tr_occ, n_tr_occ};

            auto block_g_dabi_lt_bji = g_cjkl_lt("a,j,i,n").block(g_cjkl_low_a, g_cjkl_up_a);
            auto block_g_dabi_lt_bij = g_cjkl_lt("a,i,j,n").block(g_cjkl_low_a, g_cjkl_up_a);

            T1("e,a,b,n") = block_t2_oou_lt_eb * block_g_dabi_lt_bji;
            T2("e,a,b,n") = block_t2_oou_lt_eb * block_g_dabi_lt_bij;

            /*if (this->is_df()) {

              auto n_tr_aux = Xab.range().upbound()[0];

              TArray G1;
              TArray G2;
              TArray G3;

              auto block_t2_ouu_lt_cb = t2_ouu_lt("c,b,i,n").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

              block t2_oou_lt_low_bc{b_low, 0, 0, 0};
              block t2_oou_lt_up_bc{b_up, n_tr_vir, n_tr_occ, n_tr_occ};
              auto block_t2_ouu_lt_bc = t2_ouu_lt("b,c,i,n").block(t2_oou_lt_low_bc, t2_oou_lt_up_bc);

              auto block_g_dabi_lt_ac = g_dabi_lt("e,a,c,i").block(g_dabi_low_a, g_dabi_up_a);

              block Xab_low{0, 0, a_low};
              block Xab_up{n_tr_aux, n_tr_vir, a_up};
              auto block_Xab_lt = Xab_lt("X,e,a").block(Xab_low, Xab_up);

              TArray gt1;
              gt1("b,X,n") = Xai_lt("X,c,i") * block_t2_ouu_lt_cb;
              TArray gt2;
              gt2("b,X,n") = Xai_lt("X,c,i") * block_t2_ouu_lt_bc;

              G1("e,a,b,n") = block_Xab_lt * gt1("b,X,n");
              G2("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_cb;
              G3("e,a,b,n") = block_Xab_lt * gt2("b,X,n");
              E_O2V4_2 += TA::dot(2.0*G1("e,c,a,i") - G2("e,c,a,i") - G3("e,c,a,i"),(2.0*T1("e,c,a,i") - T2("e,c,a,i")));

              TArray G;
              G("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_bc;
              E_O2V4_2 += TA::dot(G("e,c,a,i"),(-2.0*T2("e,c,a,i") + T1("e,c,a,i")));

            } else {*/

              TArray G1;
              TArray G2;
              TArray G3;

              auto block_t2_ouu_lt_cb = t2_ouu_lt("c,b,i,n").block(t2_oou_lt_low_b, t2_oou_lt_up_b);

              block t2_oou_lt_low_bc{b_low, 0, 0, 0};
              block t2_oou_lt_up_bc{b_up, n_tr_vir, n_tr_occ, n_tr_occ};
              auto block_t2_ouu_lt_bc = t2_ouu_lt("b,c,i,n").block(t2_oou_lt_low_bc, t2_oou_lt_up_bc);

              auto block_g_dabi_lt_ac = g_dabi_lt("e,a,c,i").block(g_dabi_low_a, g_dabi_up_a);

              block g_dabi_low_ca{0, 0, a_low, 0};
              block g_dabi_up_ca{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
              auto block_g_dabi_lt_ca = g_dabi_lt("e,c,a,i").block(g_dabi_low_ca, g_dabi_up_ca);

              G1("e,a,b,n") = block_g_dabi_lt_ca * block_t2_ouu_lt_cb;
              G2("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_cb;
              G3("e,a,b,n") = block_g_dabi_lt_ca * block_t2_ouu_lt_bc;
              E_O2V4_2 += TA::dot(2.0*G1("e,c,a,i") - G2("e,c,a,i") - G3("e,c,a,i"),(2.0*T1("e,c,a,i") - T2("e,c,a,i")));

              TArray G;
              G("e,a,b,n") = block_g_dabi_lt_ac * block_t2_ouu_lt_bc;
              E_O2V4_2 += TA::dot(G("e,c,a,i"),(-2.0*T2("e,c,a,i") + T1("e,c,a,i")));
            //}
          }

          {

            //if (this->is_df()) {
              /*
              auto n_tr_aux = Xab.range().upbound()[0];

              block t1_lt_low_b{b_low, 0};
              block t1_lt_up_b{b_up, n_tr_occ};
              auto block_t1_lt_b = t1_lt("b,j").block(t1_lt_low_b, t1_lt_up_b);

              auto block_t2_oou_lt_faij = t2_oou_lt("f,a,i,j").block(t2_oou_lt_low_a, t2_oou_lt_up_a);
              auto block_t2_oou_lt_faji = t2_oou_lt("f,a,j,i").block(t2_oou_lt_low_a, t2_oou_lt_up_a);
              TArray T1;
              TArray T2;
              T1("i,a,b,f") = block_t1_lt_b * block_t2_oou_lt_faij;
              T2("i,a,b,f") = block_t1_lt_b * block_t2_oou_lt_faji;

              block g_dabi_low_cb{0, 0, b_low, 0};
              block g_dabi_up_cb{n_tr_vir, n_tr_vir, b_up, n_tr_occ};
              auto block_g_dabi_lt_cb = g_dabi_lt("f,c,b,j").block(g_dabi_low_cb, g_dabi_up_cb);

              auto block_g_dabi_lt_bc = g_dabi_lt("f,b,c,j").block(g_dabi_low_b, g_dabi_up_b);

              block g_abij_lt_low_ac{a_low, 0, 0, 0};
              block g_abij_lt_up_ac{a_up, n_tr_vir, n_tr_occ, n_tr_occ};
              auto block_g_abij_lt_ac = g_abij_lt("a,c,j,i").block(g_abij_lt_low_ac, g_abij_lt_up_ac);

              block g_abij_lt_low_ca{0, a_low, 0, 0};
              block g_abij_lt_up_ca{n_tr_vir, a_up, n_tr_occ, n_tr_occ};
              auto block_g_abij_lt_ca = g_abij_lt("c,a,j,i").block(g_abij_lt_low_ca, g_abij_lt_up_ca);

              TArray gg1;
              gg1("j,X,i,Y") = Xai_lt("X,c,i") * Xai_lt("Y,c,j");
              block Xai_low{0, a_low, 0};
              block Xai_up{n_tr_aux, a_up, n_tr_occ};
              auto block_Xai_lt = Xab_lt("X,a,j").block(Xai_low, Xai_up);
              TArray ggg1;
              ggg1("a,i,Y") = gg1("j,X,i,Y") * block_Xai_lt;


              block Xab_low{0, 0, b_low};
              block Xab_up{n_tr_aux, n_tr_vir, b_up};
              auto block_Xab_lt = Xab_lt("Y,f,b").block(Xab_low, Xab_up);

              TArray G1;
              G1("i,a,b,f") = ggg1("a,i,Y") * block_Xab_lt;


              TArray G3;
              G3("i,a,b,f") = block_g_abij_lt_ac * block_g_dabi_lt_bc;
              E_O2V4_S += TA::dot((2.0*G1("i,a,b,f") - 4.0*G3("i,a,b,f")),T1("i,a,b,f"));


              TArray G2;
              G2("i,a,b,f") = block_g_abij_lt_ca * block_g_dabi_lt_bc;
              E_O2V4_S += TA::dot((-4.0*G2("i,a,b,f") - 4.0*G1("i,a,b,f")),T2("i,a,b,f"));
*/
            //} else {


              block t1_lt_low_b{b_low, 0};
              block t1_lt_up_b{b_up, n_tr_occ};
              auto block_t1_lt_b = t1_lt("b,j").block(t1_lt_low_b, t1_lt_up_b);

              auto block_t2_oou_lt_faij = t2_oou_lt("f,a,i,j").block(t2_oou_lt_low_a, t2_oou_lt_up_a);
              auto block_t2_oou_lt_faji = t2_oou_lt("f,a,j,i").block(t2_oou_lt_low_a, t2_oou_lt_up_a);
              TArray T1;
              TArray T2;
              T1("i,a,b,f") = block_t1_lt_b * block_t2_oou_lt_faij;
              T2("i,a,b,f") = block_t1_lt_b * block_t2_oou_lt_faji;

              block g_dabi_low_cb{0, 0, b_low, 0};
              block g_dabi_up_cb{n_tr_vir, n_tr_vir, b_up, n_tr_occ};
              auto block_g_dabi_lt_cb = g_dabi_lt("f,c,b,j").block(g_dabi_low_cb, g_dabi_up_cb);

              auto block_g_dabi_lt_bc = g_dabi_lt("f,b,c,j").block(g_dabi_low_b, g_dabi_up_b);

              block g_abij_lt_low_ac{a_low, 0, 0, 0};
              block g_abij_lt_up_ac{a_up, n_tr_vir, n_tr_occ, n_tr_occ};
              auto block_g_abij_lt_ac = g_abij_lt("a,c,j,i").block(g_abij_lt_low_ac, g_abij_lt_up_ac);

              block g_abij_lt_low_ca{0, a_low, 0, 0};
              block g_abij_lt_up_ca{n_tr_vir, a_up, n_tr_occ, n_tr_occ};
              auto block_g_abij_lt_ca = g_abij_lt("c,a,j,i").block(g_abij_lt_low_ca, g_abij_lt_up_ca);

              TArray G1;
              G1("i,a,b,f") = block_g_abij_lt_ac * block_g_dabi_lt_cb;

              TArray G3;
              G3("i,a,b,f") = block_g_abij_lt_ac * block_g_dabi_lt_bc;
              E_O2V4_S += TA::dot((2.0*G1("i,a,b,f") - 4.0*G3("i,a,b,f")),T1("i,a,b,f"));


              TArray G2;
              G2("i,a,b,f") = block_g_abij_lt_ca * block_g_dabi_lt_bc;
              E_O2V4_S += TA::dot((-4.0*G2("i,a,b,f") - 4.0*G1("i,a,b,f")),T2("i,a,b,f"));
            //}

            {
              block g_abij_lt_low_ab{a_low, b_low, 0, 0};
              block g_abij_lt_up_ab{a_up, b_up, n_tr_occ, n_tr_occ};
              auto block_g_abij_lt_ab = g_abij_lt("a,b,i,j").block(g_abij_lt_low_ab, g_abij_lt_up_ab);
              auto block_t2_ouu_lt_ab = t2_ouu_lt("a,b,k,n").block(g_abij_lt_low_ab, g_abij_lt_up_ab);
              TArray G;
              TArray T;
              G("a,b,c,n") = block_g_abij_lt_ab * g_cjkl_lt("c,i,j,n");
              T("a,b,c,n") = t1_lt("c,k") * block_t2_ouu_lt_ab;
              E_O3V3_S2 += 2.0 * (TA::dot((G("a,b,c,n")),(T("a,b,c,n"))));
            }
          }

        }
      }
      TA::set_default_world(global_world);
      //global_world.gop.sum(E_OV5);
      global_world.gop.sum(E_O2V4_2);
      global_world.gop.sum(E_O2V4_S);
      global_world.gop.sum(E_O3V3_S2);
      global_world.gop.sum(time_G_OV5);
      global_world.gop.sum(time_T_OV5);
      global_world.gop.sum(time_trace);


      E_O2V4_vo = 0.0;
      TA::set_default_world(this_world);
      global_iter = 0;
      for (auto a = 0; a < n_tr_vir; ++a) {

        ++ global_iter;
        if( global_iter % size != rank) continue;

        std::size_t a_low = a;
        std::size_t a_up = a + 1;

        typedef std::vector<std::size_t> block;
        {
          //blocking g_dabi over a & b
          block g_dabi_low_ab{0, a_low, 0, 0};
          block g_dabi_up_ab{n_tr_vir, a_up, n_tr_vir, n_tr_occ};
          auto block_g_dabi_lt_ab = g_dabi_lt("e,a,b,j").block(g_dabi_low_ab, g_dabi_up_ab);

          TArray G2;
          G2("e,a,i,f") = block_g_dabi_lt_ab*t2_oou_lt("f,b,j,i");
          TArray G3;
          G3("e,a,i,f") = block_g_dabi_lt_ab*t2_oou_lt("f,b,i,j");
          E_O2V4_vo += TA::dot(G2("e,a,i,f"),G2("f,a,i,e") - 4.0*G3("f,a,i,e"));
          {
            /*if (this->is_df()) {

              auto n_tr_aux = Xab.range().upbound()[0];

              block Xab_low{0, 0, a_low};
              block Xab_up{n_tr_aux, n_tr_vir, a_up};
              auto block_Xab_lt_ba = Xab_lt("X,e,a").block(Xab_low, Xab_up);

              TArray gt1;
              gt1("f,X,i") = Xai_lt("X,b,j") * t2_oou_lt("f,b,j,i");
              TArray G4;
              G4("e,a,i,f") = block_Xab_lt_ba*gt1("f,X,i");
              E_O2V4_vo += TA::dot(G4("e,a,i,f"),G4("f,a,i,e") - 4.0*G2("f,a,i,e"));
              E_O2V4_vo += TA::dot(G3("e,a,i,f"),2.0*G4("f,a,i,e") + G3("f,a,i,e"));

              TArray gt2;
              gt2("f,X,i") = Xai_lt("X,b,j") * t2_oou_lt("f,b,i,j");
              TArray G1;
              G1("e,a,i,f") = block_Xab_lt_ba*gt2("f,X,i");
              E_O2V4_vo += TA::dot(G1("e,a,i,f"),8.0*G2("f,a,i,e") + 4.0*G1("f,a,i,e") - 4.0*G4("f,a,i,e") - 4.0*G3("f,a,i,e"));

            } else {*/

              block g_dabi_low_ba{0, 0, a_low, 0};
              block g_dabi_up_ba{n_tr_vir, n_tr_vir, a_up, n_tr_occ};
              auto block_g_dabi_lt_ba = g_dabi_lt("e,b,a,j").block(g_dabi_low_ba, g_dabi_up_ba);
              TArray G4;
              G4("e,a,i,f") = block_g_dabi_lt_ba*t2_oou_lt("f,b,j,i");
              E_O2V4_vo += TA::dot(G4("e,a,i,f"),G4("f,a,i,e") - 4.0*G2("f,a,i,e"));
              E_O2V4_vo += TA::dot(G3("e,a,i,f"),2.0*G4("f,a,i,e") + G3("f,a,i,e"));
              {
                TArray G1;
                G1("e,a,i,f") = block_g_dabi_lt_ba*t2_oou_lt("f,b,i,j");
                E_O2V4_vo += TA::dot(G1("e,a,i,f"),8.0*G2("f,a,i,e") + 4.0*G1("f,a,i,e") - 4.0*G4("f,a,i,e") - 4.0*G3("f,a,i,e"));
              }

            //}
          }
        }
      }
      TA::set_default_world(global_world);
      global_world.gop.sum(E_O2V4_vo);

      //computation of the OV5 terms
      /*E_O2V4_vo = 0;
      {
        TArray G2;
        G2("e,a,i,f") = g_dabi_lt("e,a,b,j")*t2_oou_lt("f,b,j,i");
        TArray G3;
        G3("e,a,i,f") = g_dabi_lt("e,a,b,j")*t2_oou_lt("f,b,i,j");
        E_O2V4_vo = TA::dot(G2("e,a,i,f"),G2("f,a,i,e") - 4.0*G3("f,a,i,e"));
        {
          TArray G4;
          G4("e,a,i,f") = g_dabi_lt("e,b,a,j")*t2_oou_lt("f,b,j,i");
          E_O2V4_vo += TA::dot(G4("e,a,i,f"),G4("f,a,i,e") - 4.0*G2("f,a,i,e"));
          E_O2V4_vo += TA::dot(G3("e,a,i,f"),2.0*G4("f,a,i,e") + G3("f,a,i,e"));
          {
            TArray G1;
            G1("e,a,i,f") = g_dabi_lt("e,b,a,j")*t2_oou_lt("f,b,i,j");
            E_O2V4_vo += TA::dot(G1("e,a,i,f"),8.0*G2("f,a,i,e") + 4.0*G1("f,a,i,e") - 4.0*G4("f,a,i,e") - 4.0*G3("f,a,i,e"));
          }
        }
      }
      std::cout << "E_O2V4_vo2 = " << E_O2V4_vo << std::endl;*/

      //computation of the OV5 terms
      double E_O2V4_oo = 0;
      {
        TArray T1;
        T1("e,i,j,f") = t2_oou_lt("e,a,j,k")*t2_oou_lt("f,a,i,k");
        TArray T3;
        T3("e,i,j,f") = t2_oou_lt("e,a,k,j")*t2_oou_lt("f,a,k,i");
        {
          TArray T4;
          T4("e,i,j,f") = t2_oou_lt("e,a,j,k")*t2_oou_lt("f,a,k,i");
          TArray G1;
          G1("e,i,j,f") = g_dabi_lt("e,a,b,i")*g_dabi_lt("f,b,a,j");
          E_O2V4_oo = TA::dot(G1("e,i,j,f"),4.0*T1("e,i,j,f") + T3("e,i,j,f") - 4.0*T4("e,i,j,f"));
        }
        {
          TArray T2;
          T2("e,i,j,f") = t2_oou_lt("e,a,k,j")*t2_oou_lt("f,a,i,k");
          TArray G2;
          G2("e,i,j,f") = g_dabi_lt("e,a,b,i")*g_dabi_lt("f,a,b,j");
          E_O2V4_oo += TA::dot(G2("e,i,j,f"),2.0*T2("e,i,j,f") - 2.0*T3("e,i,j,f") - 2.0*T1("e,i,j,f"));
        }
      }

      //computation of the OV5 terms
      double E_OV4 = 0;
      {
        TArray T1;
        T1("e,f") = t2_oou_lt("e,a,i,j")*t2_oou_lt("f,a,i,j");
        TArray T2;
        T2("e,f") = t2_oou_lt("e,a,i,j")*t2_oou_lt("f,a,j,i");
        {
          TArray G;
          G("e,f") = g_dabi_lt("e,a,b,i")*g_dabi_lt("f,a,b,i");
          E_OV4 = TA::dot(G("e,f"),4.0*T1("e,f") - 2.0*T2("e,f"));
        }
        {
          TArray G;
          G("e,f") = g_dabi_lt("e,a,b,i")*g_dabi_lt("f,b,a,i");
          E_OV4 += TA::dot(G("e,f"),T2("e,f") - 2.0*T1("e,f"));
        }
      }

      energy_m = E_OV5 + E_O2V4_vo + E_O2V4_oo + E_OV4;

      energy_m = 6.0 * energy_m * w(m);

      time03 = mpqc::now(world, accurate_time);
      time_t2_vv += mpqc::duration_in_s(time02, time03);



      time06 = mpqc::now(world, accurate_time);

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
          E_O4V = TA::dot((G("m,n")),(4.0*T1("m,n") - 2.0*T2("m,n")));
        }
        {
          TArray G;
          TArray T;
          G("m,n") = g_cjkl_lt("a,i,j,m") * g_cjkl_lt("a,j,i,n");
          E_O4V += TA::dot((G("m,n")),(T2("m,n") - 2.0*T1("m,n")));
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
          E_O4V += TA::dot((4.0*G1("m,i,j,n") + G3("m,i,j,n") - 4.0*G4("m,i,j,n")),(T("m,i,j,n")));
        }
        {
          TArray G2;
          G2("m,i,j,n") = g_cjkl_lt("a,k,j,m") * g_cjkl_lt("a,i,k,n");
          TArray T;
          T("m,i,j,n") = t2_ouu_lt("a,b,i,m") * t2_ouu_lt("a,b,j,n");
          E_O4V += TA::dot((2.0*G2("m,i,j,n") - 2.0*G3("m,i,j,n") - 2.0*G1("m,i,j,n")),(T("m,i,j,n")));
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

          E_O4V2 = (TA::dot((G("m,a,i,n")),(8.0*T1("m,a,i,n") + T4("m,a,i,n") - 4.0*T3("m,a,i,n"))));
        }
        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,i,j,m") * t2_ouu_lt("b,a,j,n");

          E_O4V2 += (TA::dot((G("m,a,i,n")),(4.0*T1("m,a,i,n"))));
        }
        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,j,i,m") * t2_ouu_lt("b,a,j,n");

          E_O4V2 += (TA::dot((G("m,a,i,n")),(2.0*T2("m,a,i,n") + T3("m,a,i,n") - 4.0*T1("m,a,i,n"))));
        }
        {
          TArray G;
          G("m,a,i,n") = g_cjkl_lt("b,i,j,m") * t2_ouu_lt("a,b,j,n");

          E_O4V2 += (TA::dot((G("m,a,i,n")),(T2("m,a,i,n") - 4.0*T4("m,a,i,n") - 4.0*T1("m,a,i,n"))));
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

          E_O4V2 += TA::dot((G("m,a,b,n")),(4.0*T1("m,a,b,n") + T4("m,a,b,n") - 4.0*T3("m,a,b,n")));
        }
        {
          TArray G;
          G("m,a,b,n") = g_cjkl_lt("b,i,j,m") * g_cjkl_lt("a,i,j,n");

          E_O4V2 += TA::dot((G("m,a,b,n")),(2.0*T2("m,a,b,n") - 2.0*T4("m,a,b,n") - 2.0*T1("m,a,b,n")));
        }

      }
      energy_mo = E_O4V + E_O4V2;

      time07 = mpqc::now(world, accurate_time);
      time_t2_oo += mpqc::duration_in_s(time06, time07);

      time08 = mpqc::now(world, accurate_time);

      //double E_O2V4_vo = 0;
      //double E_O2V4_2 = 0;
      /*for (auto e = 0; e < n_tr_vir; ++e) {

        std::size_t e_low = e;
        std::size_t e_up = e + 1;

        typedef std::vector<std::size_t> block;

        //blocking g_dabi over e
        block g_dabi_low_e{e_low, 0, 0, 0};
        block g_dabi_up_e{e_up, n_tr_vir, n_tr_vir, n_tr_occ};

        auto block_g_dabi_lt_e = g_dabi_lt("e,b,c,i").block(g_dabi_low_e, g_dabi_up_e);

        //blocking t2 over e (required for both T1 and T2 intermediates)
        block t2_oou_lt_low_e{e_low, 0, 0, 0};
        block t2_oou_lt_up_e{e_up, n_tr_vir, n_tr_occ, n_tr_occ};

        auto block_t2_oou_lt_e_T1 = t2_oou_lt("e,b,i,j").block(t2_oou_lt_low_e, t2_oou_lt_up_e);

        TArray T1;
        TArray T2;

        T1("e,a,b,n") = block_t2_oou_lt_e_T1 * g_cjkl_lt("a,j,i,n");
        T2("e,a,b,n") = block_t2_oou_lt_e_T1 * g_cjkl_lt("a,i,j,n");

        {
          TArray G1;
          TArray G2;
          TArray G3;
          auto block_g_dabi_lt_e1 = g_dabi_lt("e,c,a,i").block(g_dabi_low_e, g_dabi_up_e);
          G1("e,a,b,n") = block_g_dabi_lt_e1 * t2_ouu_lt("c,b,i,n");
          auto block_g_dabi_lt_e2 = g_dabi_lt("e,a,c,i").block(g_dabi_low_e, g_dabi_up_e);
          G2("e,a,b,n") = block_g_dabi_lt_e2 * t2_ouu_lt("c,b,i,n");
          auto block_g_dabi_lt_e3 = g_dabi_lt("e,c,a,i").block(g_dabi_low_e, g_dabi_up_e);
          G3("e,a,b,n") = block_g_dabi_lt_e3 * t2_ouu_lt("b,c,i,n");
          E_O2V4_2 += TA::dot(2.0*G1("e,c,a,i") - G2("e,c,a,i") - G3("e,c,a,i"),(2.0*T1("e,c,a,i") - T2("e,c,a,i")));
        }
        {
          TArray G;
          auto block_g_dabi_lt_e = g_dabi_lt("e,a,c,i").block(g_dabi_low_e, g_dabi_up_e);
          G("e,a,b,n") = block_g_dabi_lt_e * t2_ouu_lt("b,c,i,n");
          E_O2V4_2 += TA::dot(G("e,c,a,i"),(-2.0*T2("e,c,a,i") + T1("e,c,a,i")));
        }
      }*/

      //computation of the O2V4 terms
      //double E_O2V4_2 = 0;
      /*{
        TArray T1;
        TArray T2;

        T1("e,a,b,n") = t2_oou_lt("e,b,i,j") * g_cjkl_lt("a,j,i,n");
        T2("e,a,b,n") = t2_oou_lt("e,b,i,j") * g_cjkl_lt("a,i,j,n");
        {
          TArray G1;
          TArray G2;
          TArray G3;
          G1("e,a,b,n") = g_dabi_lt("e,c,a,i") * t2_ouu_lt("c,b,i,n");
          G2("e,a,b,n") = g_dabi_lt("e,a,c,i") * t2_ouu_lt("c,b,i,n");
          G3("e,a,b,n") = g_dabi_lt("e,c,a,i") * t2_ouu_lt("b,c,i,n");
          E_O2V4_2 = TA::dot(2.0*G1("e,c,a,i") - G2("e,c,a,i") - G3("e,c,a,i"),(2.0*T1("e,c,a,i") - T2("e,c,a,i")));
        }
        {
          TArray G;
          G("e,a,b,n") = g_dabi_lt("e,a,c,i") * t2_ouu_lt("b,c,i,n");
          E_O2V4_2 += TA::dot(G("e,c,a,i"),(-2.0*T2("e,c,a,i") + T1("e,c,a,i")));
        }
      }*/

      double E_O3V3_2 = 0;
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

          E_O3V3_2 = TA::dot(3.0*G1("e,i,a,n"),(T1("e,i,a,n")));
          E_O3V3_2 += TA::dot((G1("e,i,a,n") + G2("e,i,a,n")),(T1("e,i,a,n") + 4.0*T2("e,i,a,n") - 2.0*T3("e,i,a,n") - 2.0*T4("e,i,a,n")));
        }
        {
          TArray G3;
          TArray G4;
          G3("e,i,a,n") = g_dabi_lt("e,a,b,j") * g_cjkl_lt("b,i,j,n");
          G4("e,i,a,n") = g_dabi_lt("e,b,a,j") * g_cjkl_lt("b,j,i,n");

          E_O3V3_2 += TA::dot((G3("e,i,a,n") + G4("e,i,a,n")),(T3("e,i,a,n") + T4("e,i,a,n") - 2.0*T1("e,i,a,n") - 2.0*T2("e,i,a,n")));

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

          E_O3V3_2 += TA::dot((G("e,i,j,n")),(4.0*T1("e,i,j,n") + T4("e,i,j,n") - 2.0*T3("e,i,j,n") - 2.0*T2("e,i,j,n")));
        }
        {
          TArray G;
          G("e,i,j,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("a,b,j,n");

          E_O3V3_2 += TA::dot((G("e,i,j,n")),(T2("e,i,j,n") + T3("e,i,j,n") - 2.0*T4("e,i,j,n") - 2.0*T1("e,i,j,n")));
        }
      }

      double E_O2V3_2 = 0;
      {
        TArray T1;
        TArray T2;

        T1("e,n") = t2_oou_lt("e,a,i,j") * g_cjkl_lt("a,i,j,n");
        T2("e,n") = t2_oou_lt("e,a,i,j") * g_cjkl_lt("a,j,i,n");
        {
          TArray G;
          G("e,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("a,b,i,n");
          E_O2V3_2 = TA::dot(G("e,n"),4.0*T1("e,n") - 2.0*T2("e,n"));
        }
        {
          TArray G;
          G("e,n") = g_dabi_lt("e,a,b,i") * t2_ouu_lt("b,a,i,n");
          E_O2V3_2 += TA::dot(G("e,n"),T2("e,n") - 2.0*T1("e,n"));
        }
      }

      time09 = mpqc::now(world, accurate_time);
      time_t2_ov += mpqc::duration_in_s(time08, time09);



      energy_muo = E_O2V4_2 + E_O2V3_2 + E_O3V3_2;


      time12 = mpqc::now(world, accurate_time);

      /*double E_O2V4_S = 0.0;
      TA::set_default_world(this_world);
      global_iter = 0;
      for (auto e = 0; e < n_tr_vir; ++e) {

        ++ global_iter;
        if( global_iter % size != rank) continue;

        std::size_t e_low = e;
        std::size_t e_up = e + 1;

        typedef std::vector<std::size_t> block;

        if (this->is_df()) {

          auto n_tr_aux = Xab.range().upbound()[0];

          //blocking g_dabi over e
          block g_dabi_low_e{e_low, 0, 0, 0};
          block g_dabi_up_e{e_up, n_tr_vir, n_tr_vir, n_tr_occ};

          //blocking t2 over e (required for both T1 and T2 intermediates)
          block t2_oou_lt_low_e{e_low, 0, 0, 0};
          block t2_oou_lt_up_e{e_up, n_tr_vir, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_e_T1 = t2_oou_lt("e,a,i,j").block(t2_oou_lt_low_e, t2_oou_lt_up_e);
          auto block_t2_oou_lt_e_T2 = t2_oou_lt("e,a,j,i").block(t2_oou_lt_low_e, t2_oou_lt_up_e);

          TArray T1;
          TArray T2;
          T1("i,a,b,e") = t1_lt("b,j") * block_t2_oou_lt_e_T1;
          T2("i,a,b,e") = t1_lt("b,j") * block_t2_oou_lt_e_T2;

          TArray gg1;
          gg1("i,j,X,Y") = Xai_lt("X,c,i") * Xai_lt("Y,c,j");
          TArray ggg1;
          ggg1("i,a,Y") = gg1("i,j,X,Y") * Xai_lt("X,a,j");

          block Xab_low{0,e_low,0};
          block Xab_up{n_tr_aux,e_up,n_tr_vir};
          auto block_Xab_lt = Xab_lt("Y,e,b").block(Xab_low, Xab_up);

          TArray G1;
          G1("i,a,b,e") = ggg1("i,a,Y") * block_Xab_lt;

          auto block_Xab_lt_c = Xab_lt("Y,e,c").block(Xab_low, Xab_up);
          TArray gg2;
          gg2("e,X,j,Y") = Xai_lt("X,c,j") * block_Xab_lt_c;
          TArray ggg2;
          ggg2("e,X,b") = gg2("e,X,j,Y") * Xai_lt("Y,b,j");

          TArray G2;
          G2("i,a,b,e") = ggg2("e,X,b") * Xai_lt("X,a,i");

          E_O2V4_S += TA::dot((-4.0*G2("i,a,b,e") - 4.0*G1("i,a,b,e")),T2("i,a,b,e"));

          {
            TArray G3;
            auto block_g_dabi_lt_e3 = g_dabi_lt("e,b,c,j").block(g_dabi_low_e, g_dabi_up_e);
            G3("i,a,b,e") = g_abij_lt("a,c,j,i") * block_g_dabi_lt_e3;
            E_O2V4_S += TA::dot((2.0*G1("i,a,b,e") - 4.0*G3("i,a,b,e")),T1("i,a,b,e"));
          }

        } else {
          //blocking g_dabi over e
          block g_dabi_low_e{e_low, 0, 0, 0};
          block g_dabi_up_e{e_up, n_tr_vir, n_tr_vir, n_tr_occ};

          auto block_g_dabi_lt_e = g_dabi_lt("e,c,b,j").block(g_dabi_low_e, g_dabi_up_e);

          //blocking t2 over e (required for both T1 and T2 intermediates)
          block t2_oou_lt_low_e{e_low, 0, 0, 0};
          block t2_oou_lt_up_e{e_up, n_tr_vir, n_tr_occ, n_tr_occ};

          auto block_t2_oou_lt_e_T1 = t2_oou_lt("e,a,i,j").block(t2_oou_lt_low_e, t2_oou_lt_up_e);
          auto block_t2_oou_lt_e_T2 = t2_oou_lt("e,a,j,i").block(t2_oou_lt_low_e, t2_oou_lt_up_e);

          TArray T1;
          TArray T2;
          T1("i,a,b,e") = t1_lt("b,j") * block_t2_oou_lt_e_T1;
          T2("i,a,b,e") = t1_lt("b,j") * block_t2_oou_lt_e_T2;

          TArray G1;
          G1("i,a,b,e") = g_abij_lt("a,c,j,i") * block_g_dabi_lt_e;

          {
            TArray G3;
            auto block_g_dabi_lt_e3 = g_dabi_lt("e,b,c,j").block(g_dabi_low_e, g_dabi_up_e);
            G3("i,a,b,e") = g_abij_lt("a,c,j,i") * block_g_dabi_lt_e3;
            E_O2V4_S += TA::dot((2.0*G1("i,a,b,e") - 4.0*G3("i,a,b,e")),T1("i,a,b,e"));
          }
          {
            TArray G2;
            auto block_g_dabi_lt_e2 = g_dabi_lt("e,b,c,j").block(g_dabi_low_e, g_dabi_up_e);
            G2("i,a,b,e") = g_abij_lt("c,a,j,i") * block_g_dabi_lt_e2;
            E_O2V4_S += TA::dot((-4.0*G2("i,a,b,e") - 4.0*G1("i,a,b,e")),T2("i,a,b,e"));
          }
        }
      }
      TA::set_default_world(global_world);
      global_world.gop.sum(E_O2V4_S);*/
      /*{
        TArray T1;
        TArray T2;
        T1("i,a,b,f") = t1_lt("b,j") * t2_oou_lt("f,a,i,j");
        T2("i,a,b,f") = t1_lt("b,j") * t2_oou_lt("f,a,j,i");

        TArray G1;
        G1("i,a,b,f") = g_abij_lt("a,c,j,i") * g_dabi_lt("f,c,b,j");

        {
          TArray G3;
          G3("i,a,b,f") = g_abij_lt("a,c,j,i") * g_dabi_lt("f,b,c,j");
          E_O2V4_S = TA::dot((2.0*G1("i,a,b,f") - 4.0*G3("i,a,b,f")),T1("i,a,b,f"));
        }
        {
          TArray G2;
          G2("i,a,b,f") = g_abij_lt("c,a,j,i") * g_dabi_lt("f,b,c,j");
          E_O2V4_S += TA::dot((-4.0*G2("i,a,b,f") - 4.0*G1("i,a,b,f")),T2("i,a,b,f"));
        }
      }*/

      double E_O2V3_S = 0.0;
      {
        {
          TArray T1;
          T1("f,i") = t1_lt("a,j") * t2_oou_lt("f,a,i,j");
          TArray G1;
          TArray G4;
          G1("f,i") = g_abij_lt("a,b,j,i") * g_dabi_lt("f,a,b,j");
          G4("f,i") = g_abij_lt("a,b,i,j") * g_dabi_lt("f,a,b,j");
          E_O2V3_S = TA::dot((8.0*G1("f,i") - 4.0*G4("f,i")),T1("f,i"));
        }
        {
          TArray T2;
          T2("f,i") = t1_lt("a,j") * t2_oou_lt("f,a,j,i");
          TArray G2;
          TArray G3;
          G2("f,i") = g_abij_lt("a,b,j,i") * g_dabi_lt("f,b,a,j");
          G3("f,i") = g_abij_lt("a,b,i,j") * g_dabi_lt("f,b,a,j");
          E_O2V3_S += TA::dot((2.0*G2("f,i") - 4.0*G3("f,i")),T2("f,i"));
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

        E_O2V3_S += TA::dot((2.0*G1("f,a") - G2("f,a")),(4.0*T1("f,a") - 2.0*T2("f,a")));
      }

      double E_O3V3_S = 0.0;
      {
        TArray T1;
        TArray T2;
        TArray G;
        T1("i,j,k,f") = t1_lt("a,k") * t2_oou_lt("f,a,i,j");
        T2("i,j,k,f") = t1_lt("a,k") * t2_oou_lt("f,a,j,i");
        G("i,j,k,f") = g_abij_lt("a,b,i,j") * g_dabi_lt("f,a,b,k");

        E_O3V3_S = TA::dot((G("i,j,k,f")),(2.0*T1("i,j,k,f") - 4.0*T2("i,j,k,f")));
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
          E_O3V3_S += TA::dot((8.0*G1("a,i,j,f") - 4.0*G3("a,i,j,f")),T1("a,i,j,f"));
        }
        {
          TArray G2;
          G2("a,i,j,f") = g_abij_lt("a,b,k,i") * t2_oou_lt("f,b,j,k");
          E_O3V3_S += TA::dot((2.0*G2("a,i,j,f") + 2.0*G3("a,i,j,f")),T2("a,i,j,f"));
        }
      }

      energy_ms = E_O2V4_S + E_O2V3_S + E_O3V3_S;

      double E_O3V2_S = 0.0;
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
          E_O3V2_S = TA::dot((8.0*G1("n,i") - 4.0*G4("n,i")),T1("n,i"));
        }
        {
          TArray G2;
          TArray G3;
          G2("n,i") = g_abij_lt("a,b,j,i") * t2_ouu_lt("b,a,j,n");
          G3("n,i") = g_abij_lt("a,b,i,j") * t2_ouu_lt("b,a,j,n");
          E_O3V2_S += TA::dot((2.0*G2("n,i") - 4.0*G3("n,i")),T2("n,i"));
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
          E_O3V2_S += TA::dot((G("n,a")),(8.0*T1("n,a") - 4.0*T2("n,a")));
        }
        {
          TArray G;
          G("n,a") = g_abij_lt("a,b,i,j") * g_cjkl_lt("b,j,i,n");
          E_O3V2_S += TA::dot((G("n,a")),(2.0*T2("n,a") - 4.0*T1("n,a")));
        }
      }

      double E_O4V2_S = 0.0;
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
          E_O4V2_S = TA::dot((8.0*G1("i,j,a,n") + 2.0*G2("i,j,a,n") - 4.0*G4("i,j,a,n")),(T1("i,j,a,n")));
        }
        {
          TArray T2;
          T2("i,j,a,n") = t1_lt("b,j") * t2_ouu_lt("b,a,i,n");
          TArray G3;
          G3("i,j,a,n") = g_abij_lt("a,b,i,k") * g_cjkl_lt("b,j,k,n");
          E_O4V2_S += TA::dot((-4.0*G3("i,j,a,n") - 4.0*G2("i,j,a,n")),(T2("i,j,a,n")));
        }
      }
      {
        TArray G;
        TArray T;
        G("i,j,k,n") = g_abij_lt("a,b,i,j") * t2_ouu_lt("a,b,k,n");
        T("i,j,k,n") = t1_lt("c,k") * g_cjkl_lt("c,j,i,n");
        E_O4V2_S += -4.0 * (TA::dot((G("i,j,k,n")),(T("i,j,k,n"))));
      }
      {
        {
          TArray G1;
          TArray T1;
          G1("i,a,b,n") = g_abij_lt("c,a,i,j") * t2_ouu_lt("b,c,j,n");
          T1("i,a,b,n") = t1_lt("b,j") * g_cjkl_lt("a,j,i,n");
          E_O3V3_S2 += 2.0 * (TA::dot((G1("i,a,b,n")),(T1("i,a,b,n"))));
        }
        {
          TArray G2;
          TArray G3;
          TArray T2;
          G2("i,a,b,n") = g_abij_lt("a,c,i,j") * t2_ouu_lt("b,c,j,n");
          G3("i,a,b,n") = g_abij_lt("a,c,i,j") * t2_ouu_lt("c,b,j,n");
          T2("i,a,b,n") = t1_lt("b,j") * g_cjkl_lt("a,i,j,n");
          E_O3V3_S2 += TA::dot((2.0*G2("i,a,b,n") - 4.0*G3("i,a,b,n")),(T2("i,a,b,n")));
        }
      }

      time13 = mpqc::now(world, accurate_time);
      time_t1 += mpqc::duration_in_s(time12, time13);

      energy_mos = E_O3V2_S + E_O4V2_S + E_O3V3_S2;

      energy_mos = 3.0 * energy_mos;
      energy_mos = energy_mos * w(m);
      triple_energyos += energy_mos;

      energy_ms = 3.0 * energy_ms;
      energy_ms = energy_ms * w(m);
      triple_energys += energy_ms;

      energy_muo = 6.0 * energy_muo * w(m);
      triple_energyuo += energy_muo;

      energy_mo = 6.0 * energy_mo * w(m);
      triple_energyo += energy_mo;

      triple_energy +=
          energy_m + energy_mo - 2.0 * energy_muo + energy_ms - energy_mos;
    }

    if (t1.world().rank() == 0) {
      std::cout << "int_transform1: " << int_transform1 << " S \n";
      std::cout << "int_transform2: " << int_transform2 << " S \n";
      std::cout << "int_transform3: " << int_transform3 << " S \n";
      std::cout << "time_t2_vv: " << time_t2_vv << " S \n";
      std::cout << "time_t2_oo: " << time_t2_oo << " S \n";
      std::cout << "time_t2_ov: " << time_t2_ov << " S \n";
      std::cout << "time_t1: " << time_t1 << " S \n";
      std::cout << "time_G_OV5: " << time_G_OV5 << " S \n";
      std::cout << "time_T_OV5: " << time_T_OV5 << " S \n";
      std::cout << "time_trace: " << time_trace << " S \n";
      std::cout << "time_gg1: " << time_gg1 << " S \n";
      std::cout << "time_ggg1: " << time_ggg1 << " S \n";
      std::cout << "time_G1: " << time_G1 << " S \n";
      std::cout << "time_gg2: " << time_gg2 << " S \n";
      std::cout << "time_ggg2: " << time_ggg2 << " S \n";
      std::cout << "time_G2: " << time_G2 << " S \n";

    }

    triple_energy = -triple_energy / (3.0 * alpha);

    return triple_energy;
  }

  void reblock() {
    auto &lcao_factory = this->lcao_factory();
    auto &world = lcao_factory.world();

    std::size_t b_occ = occ_block_size_;
    std::size_t b_vir = unocc_block_size_;

    std::size_t occ = this->trange1_engine_->get_occ();
    std::size_t vir = this->trange1_engine_->get_vir();
    std::size_t all = this->trange1_engine_->get_all();
    std::size_t n_frozen = this->trange1_engine_->get_nfrozen();

    TA::TiledRange1 old_occ = this->trange1_engine_->get_active_occ_tr1();
    TA::TiledRange1 old_vir = this->trange1_engine_->get_vir_tr1();

    // get occupied and virtual orbitals
    auto occ_space = lcao_factory.orbital_space().retrieve(OrbitalIndex(L"i"));
    auto vir_space = lcao_factory.orbital_space().retrieve(OrbitalIndex(L"a"));

    if (reblock_) {
      auto new_tr1 =
          std::make_shared<TRange1Engine>(occ, all, b_occ, b_vir, n_frozen);

      TA::TiledRange1 new_occ = new_tr1->get_active_occ_tr1();
      TA::TiledRange1 new_vir = new_tr1->get_vir_tr1();

      mpqc::detail::parallel_print_range_info(world, new_occ, "CCSD(T) Occ");
      mpqc::detail::parallel_print_range_info(world, new_vir, "CCSD(T) Vir");

      this->set_trange1_engine(new_tr1);

      TArray occ_convert =
          array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
              world, old_occ, new_occ, 1.0);

      TArray vir_convert =
          array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
              world, old_vir, new_vir, 1.0);

      auto new_occ_space = occ_space;
      new_occ_space("k,i") = occ_space("k,j") * occ_convert("j,i");

      auto new_vir_space = vir_space;
      new_vir_space("k,a") = vir_space("k,b") * vir_convert("b,a");

      lcao_factory.orbital_space().clear();
      lcao_factory.orbital_space().add(new_occ_space);
      lcao_factory.orbital_space().add(new_vir_space);

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
      auto &tr1 = this->trange1_engine();

      // occ inner
      tr_occ_inner_ =
          tr1->compute_range(tr1->get_active_occ(), inner_block_size_);

      mpqc::detail::parallel_print_range_info(world, tr_occ_inner_,
                                              "CCSD(T) OCC Inner");

      auto occ_inner_convert =
          array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
              world, old_occ, tr_occ_inner_, 1.0);

      TArray inner_occ;
      inner_occ("k,i") = occ_space("k,j") * occ_inner_convert("j,i");
      OrbitalSpace<TArray> inner_occ_space = OrbitalSpace<TArray>(
          OrbitalIndex(L"m"), OrbitalIndex(L"κ"), inner_occ);

      lcao_factory.orbital_space().remove(OrbitalIndex(L"m"));
      lcao_factory.orbital_space().add(inner_occ_space);

      // vir inner
      tr_vir_inner_ = tr1->compute_range(vir, inner_block_size_);
      mpqc::detail::parallel_print_range_info(world, tr_vir_inner_,
                                              "CCSD(T) Vir Inner");
      auto vir_inner_convert =
          array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
              world, old_vir, tr_vir_inner_, 1.0);

      TArray inner_vir;
      inner_vir("k,a") = vir_space("k,b") * vir_inner_convert("b,a");
      OrbitalSpace<TArray> inner_vir_space = OrbitalSpace<TArray>(
          OrbitalIndex(L"a'"), OrbitalIndex(L"κ"), inner_vir);

      lcao_factory.orbital_space().remove(OrbitalIndex(L"a'"));
      lcao_factory.orbital_space().add(inner_vir_space);

      utility::print_par(world,
                         "Warning!! Using m for Inner Occupied Orbitals and a' "
                         "for Inner Virtual Orbitals! \n");
    }
  }

  void reblock_inner_t2(TArray &t2_left, TArray &t2_right) {
    auto &world = this->lcao_factory().world();

    auto vir_inner_convert =
        array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
            world, t2_left.trange().data()[0], tr_vir_inner_, 1.0);

    auto occ_inner_convert =
        array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
            world, t2_right.trange().data()[3], tr_occ_inner_, 1.0);

    t2_left("d,a,i,j") = t2_left("c,a,i,j") * vir_inner_convert("c,d");

    t2_right("a,b,i,l") = t2_right("a,b,i,j") * occ_inner_convert("j,l");
  }

  const TArray get_Xab() {
    TArray result;
    TArray sqrt =
        this->lcao_factory().ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center;
    if (reblock_inner_) {
      three_center = this->lcao_factory().compute(L"(Κ|G|a' b)");
    } else {
      three_center = this->lcao_factory().compute(L"(Κ|G|a b)");
    }
    result("K,a,b") = sqrt("K,Q") * three_center("Q,a,b");
    return result;
  }

  /// get three center integral (X|ai)
  const TArray get_Xai() {
    TArray result;
    TArray sqrt =
        this->lcao_factory().ao_factory().compute(L"(Κ|G| Λ)[inv_sqr]");
    TArray three_center = this->lcao_factory().compute(L"(Κ|G|a i)");
    result("K,a,i") = sqrt("K,Q") * three_center("Q,a,i");
    return result;
  }

  /// <ai|jk>
  const TArray get_aijk() {
    std::wstring post_fix = L"";
    if (this->is_df()) {
      post_fix = L"[df]";
    }

    if (reblock_inner_) {
      return this->lcao_factory().compute(L"<a i|G|j m>" + post_fix);
    } else {
      return this->lcao_factory().compute(L"<a i|G|j k>" + post_fix);
    }
  }

  /// <ab|ci>
  const TArray get_abci() {
    std::wstring post_fix = L"";
    if (this->is_df()) {
      post_fix = L"[df]";
    }

    if (reblock_inner_) {
      return this->lcao_factory().compute(L"<a' b|G|c i>" + post_fix);
    } else {
      return this->lcao_factory().compute(L"<a b|G|c i>" + post_fix);
    }
  }

  /// <ab|ij>
  const TArray get_abij() {
    if (this->is_df()) {
      return this->lcao_factory().compute(L"<a b|G|i j>[df]");
    } else {
      return this->lcao_factory().compute(L"<a b|G|i j>");
    }
  }

 private:
  struct ReduceBase {
    typedef double result_type;
    typedef Tile argument_type;

    std::shared_ptr<Eigen::VectorXd> vec_;
    std::size_t n_occ_;
    std::size_t n_frozen_;
    std::array<std::size_t, 6> offset_;

    ReduceBase(std::shared_ptr<Eigen::VectorXd> vec, std::size_t n_occ,
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

    CCSD_T_Reduce(std::shared_ptr<Eigen::VectorXd> vec, std::size_t n_occ,
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

    CCSD_T_ReduceSymm(std::shared_ptr<Eigen::VectorXd> vec, std::size_t n_occ,
                      std::size_t n_frozen, std::array<std::size_t, 6> offset)
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

};  // class CCSD_T

#if TA_DEFAULT_POLICY == 0
extern template class CCSD_T<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class CCSD_T<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_T_H_
