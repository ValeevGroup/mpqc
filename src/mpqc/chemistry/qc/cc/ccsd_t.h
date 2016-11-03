//
// Created by Chong Peng on 8/21/15.
//

#ifndef MPQC_CCSD_T_H_H
#define MPQC_CCSD_T_H_H

#include <mpqc/chemistry/qc/cc/ccsd.h>
#include <mpqc/util/misc/print.h>

namespace mpqc {
namespace cc {

//

/**
 *  \brief CCSD_T class that compute CCSD(T) triple calculation
 *
 */

template <typename Tile, typename Policy>
class CCSD_T : public CCSD<Tile, Policy> {
 public:
  //  using Tile = TA::TensorD;
  //  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;

 private:
  bool reblock_;
  bool reblock_inner_;
  std::size_t occ_block_size_;
  std::size_t unocc_block_size_;
  std::size_t inner_block_size_;
  std::size_t increase_;
  double triples_energy_;
  TA::TiledRange1 tr_occ_inner_;
  TA::TiledRange1 tr_vir_inner_;

 public:
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
   * | reblock_inner | int | none | block size to reblock inner dimension |
   * | increase | int | 2 | number of block in virtual dimension to load at each virtual loop |
   */

  CCSD_T(const KeyVal &kv) : CCSD<Tile, Policy>(kv) {
    reblock_ = kv.exists("reblock_occ") || kv.exists("reblock_unocc");
    reblock_inner_ = kv.exists("reblock_inner");
    occ_block_size_ = kv.value<int>("reblock_occ", 8);
    unocc_block_size_ = kv.value<int>("reblock_unocc", 8);
    inner_block_size_ = kv.value<int>("reblock_inner", 128);
    increase_ = kv.value<int>("increase", 2);
  }

  virtual ~CCSD_T() {}

  double value() override {
    if (this->energy_ == 0.0) {
      auto &world = this->lcao_factory().world();

      double ccsd_corr = 0.0;

      auto time0 = mpqc::fenced_now(world);
      ccsd_corr = CCSD<Tile, Policy>::value();
      auto time1 = mpqc::fenced_now(world);
      auto duration0 = mpqc::duration_in_s(time0, time1);
      if (world.rank() == 0) {
        std::cout << "CCSD Time " << duration0 << std::endl;
      }

      time0 = mpqc::fenced_now(world);
      // clean all LCAO integral
      this->lcao_factory().registry().purge(world);

      if (reblock_) {
        reblock();
      }

      // start CCSD(T)
      if (world.rank() == 0) {
        std::cout << "\nBegining CCSD(T) " << std::endl;
      }
      TArray t1 = this->t1();
      TArray t2 = this->t2();
      triples_energy_ = compute_ccsd_t(t1, t2);
      time1 = mpqc::fenced_now(world);
      auto duration1 = mpqc::duration_in_s(time0, time1);

      if (world.rank() == 0) {
        std::cout << std::setprecision(15);
        std::cout << "(T) Energy      " << triples_energy_ << " Time " << duration1
                  << std::endl;
        //                    std::cout << "(T) Energy      " << ccsd_t_d << "
        //                    Time " << duration2 << std::endl;
        std::cout << "CCSD(T) Energy  " << triples_energy_ + ccsd_corr << std::endl;
        //                    std::cout << "CCSD(T) Energy  " << ccsd_t_d +
        //                    ccsd_corr << std::endl;
      }

      this->energy_ = ccsd_corr + triples_energy_;
    }
    return this->energy_;
  }

  void obsolete() override {
    triples_energy_ = 0.0;
    CCSD<Tile, Policy>::obsolete();
  }

 private:
  double compute_ccsd_t(TArray &t1, TArray &t2) {
    bool df = this->is_df();
    auto &world = t1.world();
    bool accurate_time = this->lcao_factory().accurate_time();

    if (df && world.rank() == 0) {
      std::cout << "Use Density Fitting Expression to avoid storing G_vovv"
                << std::endl;
    }
    // get integral
    TArray g_cjkl = get_aijk();
    TArray g_abij = get_abij();

    TArray g_dabi;
    TArray Xdb;
    TArray Xai;

    if (df) {
      Xdb = get_Xab();
      Xai = this->get_Xai();
    } else {
      g_dabi = get_abci();
    }

    // T2
    TArray t2_left = t2;
    TArray t2_right = t2;

    if (reblock_inner_) {
      reblock_inner_t2(t2_left, t2_right);
    }

    // get trange1
    auto tr_occ = this->trange1_engine_->get_occ_tr1();
    auto tr_vir = this->trange1_engine_->get_vir_tr1();

    auto n_tr_occ = this->trange1_engine_->get_occ_blocks();
    auto n_tr_vir = this->trange1_engine_->get_vir_blocks();
    auto n_tr_occ_inner = n_tr_occ;
    auto n_tr_vir_inner = n_tr_vir;
    if (reblock_inner_) {
      n_tr_occ_inner = tr_occ_inner_.tiles_range().second;
      n_tr_vir_inner = tr_vir_inner_.tiles_range().second;
    }
    std::size_t n_tr_x = 0;
    if (df) {
      n_tr_x = Xdb.trange().data().front().tiles_range().second;
    }

    double triple_energy = 0.0;

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
                 (std::pow(1024.0, 3));

    if (t1.world().rank() == 0) {
      std::cout << "Increase in the loop " << increase << std::endl;
      std::cout << "Number of blocks at each iteration " << n_blocks
                << std::endl;
      std::cout << std::setprecision(5);
      std::cout << "Size of T3 or V3 at each iteration " << mem << " GB"
                << std::endl;
    }

    double t3_time = 0.0;
    double v3_time = 0.0;
    double reduce_time = 0.0;
    double block_time = 0.0;
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

          typedef std::vector<std::size_t> block;

          time0 = mpqc::now(world, accurate_time);
          // compute t3
          TArray t3;
          // abcijk contribution
          // g^{da}_{bi}*t^{cd}_{kj} - g^{cj}_{kl}*t^{ab}_{il}
          {
            TArray block_g_dabi, block_t2_dcjk, block_g_cjkl, block_t2_abil;
            if (df) {
              // block for Xdb
              block Xdb_low{0, 0, b_low};
              block Xdb_up{n_tr_x, n_tr_vir_inner, b_up};

              // block for Xai
              block Xai_low{0, a_low, 0};
              block Xai_up{n_tr_x, a_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dabi("d,a,b,i") = Xai("X,a,i").block(Xai_low, Xai_up) *
                                        Xdb("X,d,b").block(Xdb_low, Xdb_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            } else {
              // block for g_dbai
              block g_dabi_low{0, a_low, b_low, 0};
              block g_dabi_up{n_tr_vir_inner, a_up, b_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dabi("d,a,b,i") =
                  g_dabi("d,a,b,i").block(g_dabi_low, g_dabi_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            }

            // block for t2_cdk
            block t2_dcjk_low{0, c_low, 0, 0};
            block t2_dcjk_up{n_tr_vir_inner, c_up, n_tr_occ, n_tr_occ};
            time00 = mpqc::now(world, accurate_time);
            block_t2_dcjk("d,c,j,k") =
                t2_left("d,c,j,k").block(t2_dcjk_low, t2_dcjk_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for g_cjkl
            block g_cjkl_low{c_low, 0, 0, 0};
            block g_cjkl_up{c_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_g_cjkl("c,j,k,l") =
                g_cjkl("c,j,k,l").block(g_cjkl_low, g_cjkl_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for t2_abil
            block t2_abil_low{a_low, b_low, 0, 0};
            block t2_abil_up{a_up, b_up, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_t2_abil("a,b,i,l") =
                t2_right("a,b,i,l").block(t2_abil_low, t2_abil_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            time00 = mpqc::now(world, accurate_time);
            t3("a,b,i,c,j,k") =
                block_g_dabi("d,a,b,i") * block_t2_dcjk("d,c,j,k") -
                block_t2_abil("a,b,i,l") * block_g_cjkl("c,j,k,l");
            time01 = mpqc::now(world, accurate_time);
            contraction_time1 += mpqc::duration_in_s(time00, time01);
          }

          // acbikj contribution
          // g^{da}_{ci}*t^{db}_{kj} - g^{bk}_{jl}*t^{ac}_{il}
          {
            TArray block_g_daci, block_t2_dbkj, block_g_bkjl, block_t2_acil;

            if (df) {
              // block for Xdc
              block Xdc_low{0, 0, c_low};
              block Xdc_up{n_tr_x, n_tr_vir_inner, c_up};

              // block for Xai
              block Xai_low{0, a_low, 0};
              block Xai_up{n_tr_x, a_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_daci("d,a,c,i") = Xai("X,a,i").block(Xai_low, Xai_up) *
                                        Xdb("X,d,c").block(Xdc_low, Xdc_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);

            } else {
              // block for g_daci
              block g_daci_low{0, a_low, c_low, 0};
              block g_daci_up{n_tr_vir_inner, a_up, c_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_daci("d,a,c,i") =
                  g_dabi("d,a,c,i").block(g_daci_low, g_daci_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            }

            // block for t2_bdjk
            block t2_dbkj_low{0, b_low, 0, 0};
            block t2_dbki_up{n_tr_vir_inner, b_up, n_tr_occ, n_tr_occ};

            time00 = mpqc::now(world, accurate_time);
            block_t2_dbkj("d,b,k,j") =
                t2_left("d,b,k,j").block(t2_dbkj_low, t2_dbki_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);
            // block for g_kjlb
            block g_bkjl_low{b_low, 0, 0, 0};
            block g_bkjl_up{b_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_g_bkjl("b,k,j,l") =
                g_cjkl("b,k,j,l").block(g_bkjl_low, g_bkjl_up);

            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);
            // block for t2_acil
            block t2_acil_low{a_low, c_low, 0, 0};
            block t2_acil_up{a_up, c_up, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_t2_acil("a,c,i,l") =
                t2_right("a,c,i,l").block(t2_acil_low, t2_acil_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            t3("a,c,i,b,k,j") = t3("a,b,i,c,j,k");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);
            t3("a,c,i,b,k,j") +=
                block_g_daci("d,a,c,i") * block_t2_dbkj("d,b,k,j") -
                block_t2_acil("a,c,i,l") * block_g_bkjl("b,k,j,l");
            time01 = mpqc::now(world, accurate_time);
            contraction_time2 += mpqc::duration_in_s(time00, time01);
          }

          // cabkij contribution
          // g^{dc}_{ak}*t^{db}_{ij} - g^{bi}_{jl}*t^{ca}_{kl}
          {
            TArray block_g_dcak, block_g_bijl, block_t2_dbij, block_t2_cakl;

            if (df) {
              // block for Xda
              block Xda_low{0, 0, a_low};
              block Xda_up{n_tr_x, n_tr_vir_inner, a_up};

              // block for Xck
              block Xck_low{0, c_low, 0};
              block Xck_up{n_tr_x, c_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dcak("d,c,a,k") = Xai("X,c,k").block(Xck_low, Xck_up) *
                                        Xdb("X,d,a").block(Xda_low, Xda_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);

            } else {
              // block for g_dcak
              block g_dcak_low{0, c_low, a_low, 0};
              block g_dcak_up{n_tr_vir_inner, c_up, a_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dcak("d,c,a,k") =
                  g_dabi("d,c,a,k").block(g_dcak_low, g_dcak_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            }

            // block for t2_dbij
            block t2_dbij_low{0, b_low, 0, 0};
            block t2_dbij_up{n_tr_vir_inner, b_up, n_tr_occ, n_tr_occ};
            time00 = mpqc::now(world, accurate_time);
            block_t2_dbij("d,b,i,j") =
                t2_left("d,b,i,j").block(t2_dbij_low, t2_dbij_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for g_bijl
            block g_bijl_low{b_low, 0, 0, 0};
            block g_bijl_up{b_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_g_bijl("b,i,j,l") =
                g_cjkl("b,i,j,l").block(g_bijl_low, g_bijl_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for t2_cakl
            block t2_cakl_low{c_low, a_low, 0, 0};
            block t2_cakl_up{c_up, a_up, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_t2_cakl("c,a,k,l") =
                t2_right("c,a,k,l").block(t2_cakl_low, t2_cakl_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            t3("c,a,k,b,i,j") = t3("a,c,i,b,k,j");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);

            t3("c,a,k,b,i,j") +=
                block_g_dcak("d,c,a,k") * block_t2_dbij("d,b,i,j") -
                block_t2_cakl("c,a,k,l") * block_g_bijl("b,i,j,l");
            time01 = mpqc::now(world, accurate_time);
            contraction_time3 += mpqc::duration_in_s(time00, time01);
          }

          // cbakji contribution
          // g^{dc}_{bk}*t^{da}_{ji} - g^{aj}_{il}*t^{cb}_{kl}
          {
            TArray block_g_dcbk, block_g_ajil, block_t2_daji, block_t2_cbkl;

            if (df) {
              // block for Xdb
              block Xdb_low{0, 0, b_low};
              block Xdb_up{n_tr_x, n_tr_vir_inner, b_up};

              // block for Xck
              block Xck_low{0, c_low, 0};
              block Xck_up{n_tr_x, c_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dcbk("d,c,b,k") = Xai("X,c,k").block(Xck_low, Xck_up) *
                                        Xdb("X,d,b").block(Xdb_low, Xdb_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);

            } else {
              // block for g_dcbk
              block g_dcbk_low{0, c_low, b_low, 0};
              block g_dcbk_up{n_tr_vir_inner, c_up, b_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dcbk("d,c,b,k") =
                  g_dabi("d,c,b,k").block(g_dcbk_low, g_dcbk_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            }
            // block for t2_adi
            block t2_daji_low{0, a_low, 0, 0};
            block t2_daji_up{n_tr_vir_inner, a_up, n_tr_occ, n_tr_occ};

            time00 = mpqc::now(world, accurate_time);
            block_t2_daji("d,a,j,i") =
                t2_left("d,a,j,i").block(t2_daji_low, t2_daji_up);

            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);
            // block for g_ajil
            block g_ajil_low{a_low, 0, 0, 0};
            block g_ajil_up{a_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_g_ajil("a,j,i,l") =
                g_cjkl("a,j,i,l").block(g_ajil_low, g_ajil_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for t2_cbkl
            block t2_cbkl_low{c_low, b_low, 0, 0};
            block t2_cbkl_up{c_up, b_up, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_t2_cbkl("c,b,k,l") =
                t2_right("c,b,k,l").block(t2_cbkl_low, t2_cbkl_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            t3("c,b,k,a,j,i") = t3("c,a,k,b,i,j");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);
            t3("c,b,k,a,j,i") +=
                block_g_dcbk("d,c,b,k") * block_t2_daji("d,a,j,i") -
                block_t2_cbkl("c,b,k,l") * block_g_ajil("a,j,i,l");
            time01 = mpqc::now(world, accurate_time);
            contraction_time4 += mpqc::duration_in_s(time00, time01);
          }

          // bcajki contribution
          // g^{db}_{cj}*t^{da}_{ki} - g^{ak}_{il}*t^{bc}_{jl}
          {
            TArray block_g_dbcj, block_g_akil, block_t2_daki, block_t2_bcjl;

            if (df) {
              // block for Xdc
              block Xdc_low{0, 0, c_low};
              block Xdc_up{n_tr_x, n_tr_vir_inner, c_up};

              // block for Xbj
              block Xbj_low{0, b_low, 0};
              block Xbj_up{n_tr_x, b_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dbcj("d,b,c,j") = Xai("X,b,j").block(Xbj_low, Xbj_up) *
                                        Xdb("X,d,c").block(Xdc_low, Xdc_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            } else {
              // block for g_djcb
              block g_dbcj_low{0, b_low, c_low, 0};
              block g_dbcj_up{n_tr_vir_inner, b_up, c_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dbcj("d,j,c,b") =
                  g_dabi("d,j,c,b").block(g_dbcj_low, g_dbcj_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            }

            // block for t2_daki
            block t2_daki_low{0, a_low, 0, 0};
            block t2_daki_up{n_tr_vir_inner, a_up, n_tr_occ, n_tr_occ};
            time00 = mpqc::now(world, accurate_time);
            block_t2_daki("d,a,k,i") =
                t2_left("d,a,k,i").block(t2_daki_low, t2_daki_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for g_akil
            block g_akil_low{a_low, 0, 0, 0};
            block g_akil_up{a_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_g_akil("a,k,i,l") =
                g_cjkl("a,k,i,l").block(g_akil_low, g_akil_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);
            // block for t2_bcjl
            block t2_bcjl_low{b_low, c_low, 0, 0};
            block t2_bcjl_up{b_up, c_up, n_tr_occ, n_tr_occ_inner};
            time00 = mpqc::now(world, accurate_time);
            block_t2_bcjl("b,c,j,l") =
                t2_right("b,c,j,l").block(t2_bcjl_low, t2_bcjl_up);

            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            t3("b,c,j,a,k,i") = t3("c,b,k,a,j,i");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);
            t3("b,c,j,a,k,i") +=
                block_g_dbcj("d,b,c,j") * block_t2_daki("d,a,k,i") -
                block_t2_bcjl("b,c,j,l") * block_g_akil("a,k,i,l");
            time01 = mpqc::now(world, accurate_time);
            contraction_time5 += mpqc::duration_in_s(time00, time01);
          }

          // bacjik contribution
          // g^{db}_{aj}*t^{dc}_{ik} - g^{ci}_{kl}*t^{ba}_{jl}
          {
            TArray block_g_dbaj, block_t2_dcik, block_g_cikl, block_t2_bajl;

            if (df) {
              // block for Xda
              block Xda_low{0, 0, a_low};
              block Xda_up{n_tr_x, n_tr_vir_inner, a_up};

              // block for Xbj
              block Xbj_low{0, b_low, 0};
              block Xbj_up{n_tr_x, b_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dbaj("d,b,a,j") = Xai("X,b,j").block(Xbj_low, Xbj_up) *
                                        Xdb("X,d,a").block(Xda_low, Xda_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);

            } else {
              // block for g_dbaj
              block g_dbaj_low{0, b_low, a_low, 0};
              block g_dbaj_up{n_tr_vir_inner, b_up, a_up, n_tr_occ};

              time00 = mpqc::now(world, accurate_time);
              block_g_dbaj("d,b,a,j") =
                  g_dabi("d,b,a,j").block(g_dbaj_low, g_dbaj_up);
              time01 = mpqc::now(world, accurate_time);
              block_time += mpqc::duration_in_s(time00, time01);
            }

            // block for t2_cdki
            block t2_dcik_low{0, c_low, 0, 0};
            block t2_dcik_up{n_tr_vir_inner, c_up, n_tr_occ, n_tr_occ};

            time00 = mpqc::now(world, accurate_time);
            block_t2_dcik("d,c,i,k") =
                t2_left("d,c,i,k").block(t2_dcik_low, t2_dcik_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for g_iklc
            block g_cikl_low{c_low, 0, 0, 0};
            block g_cikl_up{c_up, n_tr_occ, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_g_cikl("c,i,k,l") =
                g_cjkl("c,i,k,l").block(g_cikl_low, g_cikl_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            // block for t2_bajl
            block t2_bajl_low{b_low, a_low, 0, 0};
            block t2_bajl_up{b_up, a_up, n_tr_occ, n_tr_occ_inner};

            time00 = mpqc::now(world, accurate_time);
            block_t2_bajl("b,a,j,l") =
                t2_right("b,a,j,l").block(t2_bajl_low, t2_bajl_up);
            time01 = mpqc::now(world, accurate_time);
            block_time += mpqc::duration_in_s(time00, time01);

            t3("b,a,j,c,i,k") = t3("b,c,j,a,k,i");
            time00 = mpqc::now(world, accurate_time);
            t3_permute_time += mpqc::duration_in_s(time01, time00);
            t3("b,a,j,c,i,k") +=
                block_g_dbaj("d,b,a,j") * block_t2_dcik("d,c,i,k") -
                block_t2_bajl("b,a,j,l") * block_g_cikl("c,i,k,l");
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
            // block for g_bcjk
            block g_bcjk_low{b_low, c_low, 0, 0};
            block g_bcjk_up{b_up, c_up, n_tr_occ, n_tr_occ};

            // block for t1_ai
            block t1_ai_low{a_low, 0};
            block t1_ai_up{a_up, n_tr_occ};

            time00 = mpqc::now(world, accurate_time);
            v3("b,c,j,k,a,i") = g_abij("b,c,j,k").block(g_bcjk_low, g_bcjk_up) *
                                t1("a,i").block(t1_ai_low, t1_ai_up);
            time01 = mpqc::now(world, accurate_time);
            v3_contraction_time += mpqc::duration_in_s(time00, time01);
          }

          // acbikj contribution
          // g^{ac}_{ik}*t^{b}_{j}
          {
            // block for g_acik
            block g_acik_low{a_low, c_low, 0, 0};
            block g_acik_up{a_up, c_up, n_tr_occ, n_tr_occ};

            // block for t1_bj
            block t1_bj_low{b_low, 0};
            block t1_bj_up{b_up, n_tr_occ};

            time00 = mpqc::now(world, accurate_time);
            v3("a,c,i,k,b,j") = v3("b,c,j,k,a,i");
            time01 = mpqc::now(world, accurate_time);
            v3_permute_time += mpqc::duration_in_s(time00, time01);

            v3("a,c,i,k,b,j") +=
                g_abij("a,c,i,k").block(g_acik_low, g_acik_up) *
                t1("b,j").block(t1_bj_low, t1_bj_up);

            time00 = mpqc::now(world, accurate_time);
            v3_contraction_time += mpqc::duration_in_s(time01, time00);
          }

          // abcijk contribution
          // g^{ab}_{ij}*t^{c}_{k}
          {
            // block for g_abij
            block g_abij_low{a_low, b_low, 0, 0};
            block g_abij_up{a_up, b_up, n_tr_occ, n_tr_occ};

            // block for t1_ck
            block t1_ck_low{c_low, 0};
            block t1_ck_up{c_up, n_tr_occ};

            time00 = mpqc::now(world, accurate_time);
            v3("a,b,i,j,c,k") = v3("a,c,i,k,b,j");
            time01 = mpqc::now(world, accurate_time);
            v3_permute_time += mpqc::duration_in_s(time00, time01);
            v3("a,b,i,j,c,k") +=
                g_abij("a,b,i,j").block(g_abij_low, g_abij_up) *
                t1("c,k").block(t1_ck_low, t1_ck_up);
            time00 = mpqc::now(world, accurate_time);
            v3_contraction_time += mpqc::duration_in_s(time01, time00);
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
        print_progress(a, a + increase, n_tr_vir);
      }
      a += a_increase;
    }

    if (t1.world().rank() == 0) {
      std::cout << "Total Blocks Computed  " << n_blocks_computed;
      std::cout << " from " << std::pow(n_tr_occ, 3) * std::pow(n_tr_vir, 3)
                << std::endl;

      std::cout << "T3 Total Time: " << t3_time << " S \n";
      std::cout << "T3 Block Time: " << block_time << " S \n";
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

    TA::TiledRange1 old_occ = this->trange1_engine_->get_occ_tr1();
    TA::TiledRange1 old_vir = this->trange1_engine_->get_vir_tr1();

    auto new_tr1 =
        std::make_shared<TRange1Engine>(occ, all, b_occ, b_vir, n_frozen);

    TA::TiledRange1 new_occ = new_tr1->get_occ_tr1();
    TA::TiledRange1 new_vir = new_tr1->get_vir_tr1();

    detail::parallel_print_range_info(world, new_occ, "CCSD(T) Occ");
    detail::parallel_print_range_info(world, new_vir, "CCSD(T) Vir");

    this->set_trange1_engine(new_tr1);

    TArray occ_convert = array_ops::create_diagonal_array_from_eigen<Tile>(
        world, old_occ, new_occ, 1.0);

    TArray vir_convert = array_ops::create_diagonal_array_from_eigen<Tile>(
        world, old_vir, new_vir, 1.0);

    // get occupied and virtual orbitals
    auto occ_space = lcao_factory.orbital_space().retrieve(OrbitalIndex(L"i"));
    auto vir_space = lcao_factory.orbital_space().retrieve(OrbitalIndex(L"a"));

    auto new_occ_space = occ_space;
    new_occ_space("k,i") = occ_space("k,j") * occ_convert("j,i");

    auto new_vir_space = vir_space;
    new_vir_space("k,a") = vir_space("k,b") * vir_convert("b,a");

    lcao_factory.orbital_space().clear();
    lcao_factory.orbital_space().add(new_occ_space);
    lcao_factory.orbital_space().add(new_vir_space);

    if (reblock_inner_) {
      // occ inner
      tr_occ_inner_ =
          new_tr1->compute_range(new_tr1->get_active_occ(), inner_block_size_);

      detail::parallel_print_range_info(world, tr_occ_inner_,
                                         "CCSD(T) OCC Inner");

      auto occ_inner_convert =
          array_ops::create_diagonal_array_from_eigen<Tile>(world, old_occ,
                                                            tr_occ_inner_, 1.0);

      TArray inner_occ;
      inner_occ("k,i") = occ_space("k,j") * occ_inner_convert("j,i");
      OrbitalSpace<TArray> inner_occ_space = OrbitalSpace<TArray>(
          OrbitalIndex(L"m"), OrbitalIndex(L"κ"), inner_occ);

      lcao_factory.orbital_space().add(inner_occ_space);

      // vir inner
      tr_vir_inner_ = new_tr1->compute_range(vir, inner_block_size_);
      detail::parallel_print_range_info(world, tr_vir_inner_,
                                         "CCSD(T) Vir Inner");
      auto vir_inner_convert =
          array_ops::create_diagonal_array_from_eigen<Tile>(world, old_vir,
                                                            tr_vir_inner_, 1.0);

      TArray inner_vir;
      inner_vir("k,a") = vir_space("k,b") * vir_inner_convert("b,a");
      OrbitalSpace<TArray> inner_vir_space = OrbitalSpace<TArray>(
          OrbitalIndex(L"a'"), OrbitalIndex(L"κ"), inner_vir);
      lcao_factory.orbital_space().add(inner_vir_space);

      utility::print_par(world,
                         "Warning!! Using m for Inner Occupied Orbitals and a' "
                         "for Inner Virtual Orbitals! \n");
    }

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

  void reblock_inner_t2(TArray &t2_left, TArray &t2_right) {
    auto &world = this->lcao_factory().world();

    auto vir_inner_convert = array_ops::create_diagonal_array_from_eigen<Tile>(
        world, t2_left.trange().data()[0], tr_vir_inner_, 1.0);

    auto occ_inner_convert = array_ops::create_diagonal_array_from_eigen<Tile>(
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

}  // namespace cc
}  // namespace mpqc

#endif  // MPQC_CCSD_T_H_H
