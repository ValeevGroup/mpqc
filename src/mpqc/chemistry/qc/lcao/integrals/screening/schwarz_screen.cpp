#include "mpqc/chemistry/qc/lcao/integrals/screening/schwarz_screen.h"
#include <mpqc/util/misc/exenv.h>

#include "mpqc/math/groups/petite_list.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

int64_t Qmatrix::f2s(Qmatrix::f2s_map const &map, int64_t ind) const {
  auto it = map.find(ind);
  // If this hits the issue was likely an index that was not the
  // first function in the shell.
  assert(it != map.end());
  return it->second;
}

double Qmatrix::operator()() const { return max_elem_in_Q_; }

double Qmatrix::operator()(int64_t a) const {
  return max_elem_in_row_[f2s(f2s_maps_[0], a)];
}

double Qmatrix::max_in_row(int64_t a) const { return max_elem_in_row_[a]; }

double Qmatrix::operator()(int64_t a, int64_t b) const {
  return Q_(f2s(f2s_maps_[0], a), f2s(f2s_maps_[1], b));
}

bool Qmatrix::is_aux_Q() const { return is_aux_; }

SchwarzScreen::SchwarzScreen(std::shared_ptr<Qmatrix> Qab,
                             std::shared_ptr<Qmatrix> Qcd, double thresh)
    : Screener(),
      Qab_(std::move(Qab)),
      Qcd_(std::move(Qcd)),
      thresh_(thresh),
      thresh2_(thresh_ * thresh_) {}

double SchwarzScreen::skip_threshold() const { return thresh_; }

bool SchwarzScreen::skip(int64_t a) { return skip_(a); }
bool SchwarzScreen::skip(int64_t a) const { return skip_(a); }

bool SchwarzScreen::skip(int64_t a, int64_t b) { return skip_(a, b); }
bool SchwarzScreen::skip(int64_t a, int64_t b) const { return skip_(a, b); }

bool SchwarzScreen::skip(int64_t a, int64_t b, int64_t c) {
  return skip_(a, b, c);
}
bool SchwarzScreen::skip(int64_t a, int64_t b, int64_t c) const {
  return skip_(a, b, c);
}

bool SchwarzScreen::skip(int64_t a, int64_t b, int64_t c, int64_t d) {
  return skip_(a, b, c, d);
}
bool SchwarzScreen::skip(int64_t a, int64_t b, int64_t c, int64_t d) const {
  return skip_(a, b, c, d);
}

boost::optional<double> SchwarzScreen::estimate(int64_t a) const {
  return Qab()(a) * Qcd()();
}

boost::optional<double> SchwarzScreen::estimate(int64_t a, int64_t b) const {
  if (Qab_->is_aux_Q()) {
    return Qab()(a) * Qcd()(b);
  } else {
    return Qab()(a, b) * Qcd()();
  }
}

boost::optional<double> SchwarzScreen::estimate(int64_t a, int64_t b,
                                                int64_t c) const {
  if (Qab_->is_aux_Q()) {
    return Qab()(a) * Qcd()(b, c);
  } else {
    return Qab()(a, b) * Qcd()(c);
  }
}

boost::optional<double> SchwarzScreen::estimate(int64_t a, int64_t b, int64_t c,
                                                int64_t d) const {
  assert(!Qab_->is_aux_Q());
  return Qab()(a, b) * Qcd()(c, d);
}

namespace detail {
class QTileEstimator {
 public:
  void operator()(std::array<unsigned long, 3> sh,
                  std::array<unsigned long, 3> nsh, float *out_value,
                  float scaling_factor, Qmatrix const *Qa,
                  Qmatrix const *Qbc) const {
    auto const &Ma = Qa->Q();
    auto const &Mbc = Qbc->Q();

    const auto sh0 = sh[0];
    const auto sh1 = sh[1];
    const auto sh2 = sh[2];
    const auto stop0 = sh0 + nsh[0];
    const auto stop1 = sh1 + nsh[1];
    const auto stop2 = sh2 + nsh[2];

    float bc_sum = 0.0;
    for (auto b = sh1; b < stop1; ++b) {
      for (auto c = sh2; c < stop2; ++c) {
        bc_sum += Mbc(b, c);
      }
    }

    float norm = 0.0;
    for (auto a = sh0; a < stop0; ++a) {
      norm += Ma(a) * bc_sum;
    }

    *out_value = scaling_factor * std::sqrt(norm);
  }

  void operator()(std::array<unsigned long, 4> sh,
                  std::array<unsigned long, 4> nsh, float *out_value,
                  float scaling_factor, Qmatrix const *Qab,
                  Qmatrix const *Qcd) const {
    auto const &Mab = Qab->Q();
    auto const &Mcd = Qcd->Q();

    const auto sh0 = sh[0];
    const auto sh1 = sh[1];
    const auto sh2 = sh[2];
    const auto sh3 = sh[3];
    const auto stop0 = sh0 + nsh[0];
    const auto stop1 = sh1 + nsh[1];
    const auto stop2 = sh2 + nsh[2];
    const auto stop3 = sh3 + nsh[3];

    // Precompute the CD sum
    float cd_sum = 0.0;
    for (auto c = sh2; c < stop2; ++c) {
      for (auto d = sh3; d < stop3; ++d) {
        cd_sum += Mcd(c, d);
      }
    }

    float norm = 0.0;
    for (auto a = sh0; a < stop0; ++a) {
      for (auto b = sh1; b < stop1; ++b) {
        norm += Mab(a, b) * cd_sum;
      }
    }

    *out_value = scaling_factor * std::sqrt(norm);
  }
};

}  // namespace detail

TA::Tensor<float> SchwarzScreen::norm_estimate(
    madness::World &world, std::vector<gaussian::Basis> const &bs_array,
    TA::Pmap const &pmap, const math::PetiteList &plist, bool replicate) const {
  const auto ndims = bs_array.size();
  auto trange = gaussian::detail::create_trange(bs_array);
  auto const &tile_range = trange.tiles_range();
  auto norms = TA::Tensor<float>(trange.tiles_range(), 0.0);

  if (ndims == 3) {
    auto const &bs0 = bs_array[0];
    auto const &bs1 = bs_array[1];
    auto const &bs2 = bs_array[2];

    // Get shells for each cluster
    auto const &cs0 = bs0.cluster_shells();
    auto const &cs1 = bs1.cluster_shells();
    auto const &cs2 = bs2.cluster_shells();

    const auto csize0 = cs0.size();
    const auto csize1 = cs1.size();
    const auto csize2 = cs2.size();

    unsigned long sh0 = 0;
    for (auto c0 = 0ul; c0 < csize0; ++c0) {
      unsigned long sh1 = 0;
      unsigned long nsh0 = cs0[c0].size();

      for (auto c1 = 0ul; c1 < csize1; ++c1) {
        unsigned long sh2 = 0;
        unsigned long nsh1 = cs1[c1].size();

        for (auto c2 = 0ul; c2 < csize2; ++c2) {
          unsigned long nsh2 = cs2[c2].size();

          const auto ord =
              tile_range.ordinal(std::array<decltype(c0), 3>{c0, c1, c2});
          if (pmap.is_local(ord) && plist.is_canonical(c0, c1, c2)) {
            const float multiplicity =
                static_cast<float>(plist.multiplicity(c0, c1, c2));

            detail::QTileEstimator estimator;
            std::array<unsigned long, 3> sh{sh0, sh1, sh2};
            std::array<unsigned long, 3> nsh{nsh0, nsh1, nsh2};
            auto task_f = [=](float *out_value) {
              estimator(sh, nsh, out_value, multiplicity, Qab_.get(),
                        Qcd_.get());
            };
            world.taskq.add(task_f, &norms[ord]);
          }

          sh2 += nsh2;
        }
        sh1 += nsh1;
      }
      sh0 += nsh0;
    }  // End estimate
  } else if (ndims == 4) {
    auto const &bs0 = bs_array[0];
    auto const &bs1 = bs_array[1];
    auto const &bs2 = bs_array[2];
    auto const &bs3 = bs_array[3];

    // Get shells for each cluster
    auto const &cs0 = bs0.cluster_shells();
    auto const &cs1 = bs1.cluster_shells();
    auto const &cs2 = bs2.cluster_shells();
    auto const &cs3 = bs3.cluster_shells();

    const auto csize0 = cs0.size();
    const auto csize1 = cs1.size();
    const auto csize2 = cs2.size();
    const auto csize3 = cs3.size();

    unsigned long sh0 = 0;
    for (auto c0 = 0ul; c0 < csize0; ++c0) {
      unsigned long sh1 = 0;
      unsigned long nsh0 = cs0[c0].size();

      for (auto c1 = 0ul; c1 < csize1; ++c1) {
        unsigned long sh2 = 0;
        unsigned long nsh1 = cs1[c1].size();

        for (auto c2 = 0ul; c2 < csize2; ++c2) {
          unsigned long sh3 = 0;
          unsigned long nsh2 = cs2[c2].size();

          for (auto c3 = 0ul; c3 < csize3; ++c3) {
            unsigned long nsh3 = cs3[c3].size();

            const auto ord =
                tile_range.ordinal(std::array<decltype(c0), 4>{c0, c1, c2, c3});
            if (pmap.is_local(ord) && plist.is_canonical(c0, c1, c2, c3)) {
              const float multiplicity =
                  static_cast<float>(plist.multiplicity(c0, c1, c2, c3));

              detail::QTileEstimator estimator;
              std::array<unsigned long, 4> sh{sh0, sh1, sh2, sh3};
              std::array<unsigned long, 4> nsh{nsh0, nsh1, nsh2, nsh3};
              auto task_f = [=](float *out_value) {
                estimator(sh, nsh, out_value, multiplicity, Qab_.get(),
                          Qcd_.get());
              };

              // task_f(&norms[ord]);
              world.taskq.add(task_f, &norms[ord]);
            }

            sh3 += nsh3;
          }
          sh2 += nsh2;
        }
        sh1 += nsh1;
      }
      sh0 += nsh0;
    }  // End estimate
  } else {
    norms = Screener::norm_estimate(world, bs_array, pmap, plist);
  }
  world.gop.fence();

  if (replicate) {  // construct the sum
    world.gop.sum(norms.data(), norms.size());
  }
  world.gop.fence();

  return norms;
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
