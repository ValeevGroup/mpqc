//
// Created by Chong Peng on 8/21/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_T_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_T_H_

#include "mpqc/chemistry/qc/cc/ccsd.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/misc/print.h"

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
   */
  // clang-format on

  CCSD_T(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    reblock_ = kv.exists("reblock_occ") || kv.exists("reblock_unocc");

    occ_block_size_ = kv.value<int>("reblock_occ", 8);
    unocc_block_size_ = kv.value<int>("reblock_unocc", 8);
    replicate_ = kv.value<bool>("replicate", false);

    // default value is size of total number of orbitals, which makes it 1 block
    std::size_t n_obs = this->wfn_world()->basis_registry()->retrieve(OrbitalIndex(L"κ")).nfunctions();
    inner_block_size_ = kv.value<int>("reblock_inner", n_obs);
    reblock_inner_ = (inner_block_size_ == 0) ? false : true;

    increase_ = kv.value<int>("increase", 2);
    approach_ = kv.value<std::string>("approach", "coarse");
    if (approach_ != "coarse" && approach_ != "fine" &&
        approach_ != "straight") {
      throw InputError("Invalid (T) approach! \n", __FILE__, __LINE__, "approach");
    }
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

    // reblock occ and unocc space
    if (reblock_ || reblock_inner_) {
      reblock();
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
    }
    auto time1 = mpqc::fenced_now(world);
    auto duration1 = mpqc::duration_in_s(time0, time1);

    ExEnv::out0() << "(T) Energy: " << triples_energy_ << " Time: " << duration1
                  << " S \n";
  }

  void evaluate(Energy* result) override {
    if (!this->computed()) {
      auto &world = this->lcao_factory().world();

      CCSD<Tile, Policy>::evaluate(result);
      double ccsd_energy = this->get_value(result).derivs(0)[0];

      auto time0 = mpqc::fenced_now(world);
      // compute
      compute_ccsd_t();

      this->computed_ = true;
      this->set_value(result,ccsd_energy + triples_energy_);
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
    if (reblock_inner_) {
      n_tr_occ_inner = tr_occ_inner_.tiles_range().second;
      n_tr_vir_inner = tr_vir_inner_.tiles_range().second;
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
              ((block_g_dabi * block_t2_dcjk).set_world(this_world) -
               (block_t2_abil * block_g_cjkl).set_world(this_world))
                  .set_world(this_world);
        } else {
          t3("a,b,i,c,j,k") =
              ((block_g_dabi * block_t2_dcjk).set_world(this_world) -
               (block_t2_abil * block_g_cjkl).set_world(this_world))
                  .set_world(this_world);
        }
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
