//
// Created by Chong Peng on 6/24/15.
//

#ifndef MPQC_MP2_H
#define MPQC_MP2_H

#include "../../../../../common/namespaces.h"
#include "../../../../../include/tiledarray.h"
#include "../../../../../utility/parallel_print.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>
//#include <mpqc/chemistry/qc/f12/mp2f12.h>
#include <mpqc/chemistry/qc/scf/mo_build.h>

using namespace mpqc;

namespace mpqc {
namespace mbpt {


template <typename Tile, typename Policy>
class MP2 {

//  friend class f12::MP2F12<Tile>;
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralType = integrals::MolecularIntegral<Tile, Policy>;

  MP2() = default;

  /// constructor using MO Integral with orbitals computed
  MP2(MolecularIntegralType &mo_int,
      std::shared_ptr<Eigen::VectorXd> orbital_energy,
      const std::shared_ptr<TRange1Engine> tre)
      : mo_int_(mo_int),
        orbital_energy_(orbital_energy),
        trange1_engine_(tre) {}

  /// constructfor using MO Integral without orbitals computed, call init before compute
  MP2(MolecularIntegralType &mo_int)
      : mo_int_(mo_int) {}

  MolecularIntegralType &mo_integral() const { return mo_int_; }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return trange1_engine_;
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return orbital_energy_;
  }

  /// compute function to mp2 energy
  virtual double compute(const rapidjson::Document &in) {
    // initialize
    init(in);

    std::string method =
        in.HasMember("Method") ? in["Method"].GetString() : "df";

    double mp2_energy = 0.0;

    if (method == "four center") {
      mp2_energy = compute_four_center();
    } else if (method == "df") {
      mp2_energy = compute_df();
    } else {
      throw std::runtime_error("Wrong MP2 Method");
    }

    return mp2_energy;
  }

  /// initialize orbitals
  virtual void init(const rapidjson::Document &in) {
    if(orbital_energy_== nullptr || trange1_engine_ == nullptr) {
      auto mol = mo_int_.atomic_integral().molecule();
      int occ = mol.occupation(0) / 2;
      Eigen::VectorXd orbital_energy;
      trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(mo_int_, orbital_energy, in, mol, occ);
      orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }

 protected:


  double compute_df() {
    auto g_ijab = mo_int_.compute(L"<i j|G|a b>[df]");
    // compute mp2 energy
    double energy_mp2 =
        (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
            .reduce(
                Mp2Energy(orbital_energy_,
                          trange1_engine_->get_occ(),
                          trange1_engine_->get_nfrozen())
            );

    if (g_ijab.get_world().rank() == 0) {
      std::cout << "MP2 Energy With DF: " << energy_mp2 << std::endl;
    }

    return energy_mp2;
  }

  double compute_four_center() {
    auto g_ijab = mo_int_.compute(L"<i j|G|a b>");
    // compute mp2 energy
    double energy_mp2 =
        (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
            .reduce(
                Mp2Energy(orbital_energy_, trange1_engine_->get_occ(), trange1_engine_->get_nfrozen()));

    if (g_ijab.get_world().rank() == 0) {
      std::cout << "MP2 Energy  " << energy_mp2 << std::endl;
    }

    return energy_mp2;
  }

 private:


  struct Mp2Energy {
    using result_type = double;
    using argument_type = Tile;

    std::shared_ptr<Eig::VectorXd> vec_;
    std::size_t n_occ_;
    std::size_t n_frozen_;

    Mp2Energy(std::shared_ptr<Eig::VectorXd> vec, std::size_t n_occ, std::size_t n_frozen )
        : vec_(std::move(vec)), n_occ_(n_occ), n_frozen_(n_frozen) {}

    Mp2Energy(Mp2Energy const &) = default;

    result_type operator()() const { return 0.0; }

    result_type operator()(result_type const &t) const { return t; }

    void operator()(result_type &me, result_type const &other) const {
      me += other;
    }

    void operator()(result_type &me, argument_type const &tile) const {
      auto const &range = tile.range();
      auto const &vec = *vec_;
      auto const st = range.lobound_data();
      auto const fn = range.upbound_data();
      auto tile_idx = 0;

      auto sti = st[0];
      auto fni = fn[0];
      auto stj = st[1];
      auto fnj = fn[1];
      auto sta = st[2];
      auto fna = fn[2];
      auto stb = st[3];
      auto fnb = fn[3];

      for (auto i = sti; i < fni; ++i) {
        const auto e_i = vec[i + n_frozen_];
        for (auto j = stj; j < fnj; ++j) {
          const auto e_ij = e_i + vec[j + n_frozen_];
          for (auto a = sta; a < fna; ++a) {
            const auto e_ija = e_ij - vec[a + n_occ_];
            for (auto b = stb; b < fnb; ++b, ++tile_idx) {
              const auto e_iajb = e_ija - vec[b + n_occ_];
              me += 1.0 / (e_iajb)*tile.data()[tile_idx];
            }
          }
        }
      }
    }
  };

protected:
  MolecularIntegralType &mo_int_;
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;
  std::shared_ptr<mpqc::TRange1Engine> trange1_engine_;
};

}  // end of namespace mbpt
}  // end of namespace mpqc
#endif  // MPQC_MP2_H
