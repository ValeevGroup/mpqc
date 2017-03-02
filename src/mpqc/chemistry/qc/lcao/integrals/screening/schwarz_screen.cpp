#include "mpqc/chemistry/qc/lcao/integrals/screening/schwarz_screen.h"
#include <mpqc/util/misc/exenv.h>

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

TA::Tensor<float> SchwarzScreen::norm_estimate(
    madness::World &world, std::vector<gaussian::Basis> const &bs_array,
    const math::PetiteList& plist) const {
  const auto ndims = bs_array.size();
  auto trange = gaussian::detail::create_trange(bs_array);
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

    auto sh0 = 0;
    for (auto c0 = 0ul; c0 < csize0; ++c0) {
      auto sh1 = 0;
      const auto nsh0 = cs0[c0].size();

      for (auto c1 = 0ul; c1 < csize1; ++c1) {
        auto sh2 = 0;
        const auto nsh1 = cs1[c1].size();

        for (auto c2 = 0ul; c2 < csize2; ++c2) {
          const auto nsh2 = cs2[c2].size();

          auto task_f = [=](float *out, int64_t scaling_factor) {
            auto &Qab = Qab_->Q();
            auto &Qcd = Qcd_->Q();

            float norm = 0.0;

            for (auto a = sh0; a < nsh0 + sh0; ++a) {
              const auto Qa = Qab(a);

              for (auto b = sh1; b < nsh1 + sh1; ++b) {
                for (auto c = sh2; c < nsh2 + sh2; ++c) {
                  const auto val = a * Qcd(b, c);
                  norm += val;
                }
              }
            }

            *out = scaling_factor * std::sqrt(norm);
          };

          if (plist.is_canonical(c0,c1,c2)) {
            const float multiplicity = static_cast<float>(plist.multiplicity(c0,c1,c2));
            world.taskq.add(task_f, &norms(c0, c1, c2), multiplicity);
          }

          sh2 += nsh2;
        }
        sh1 += nsh1;
      }
      sh0 += nsh0;
    } // End estimate
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

    auto sh0 = 0;
    for (auto c0 = 0ul; c0 < csize0; ++c0) {
      auto sh1 = 0;
      const auto nsh0 = cs0[c0].size();

      for (auto c1 = 0ul; c1 < csize1; ++c1) {
        auto sh2 = 0;
        const auto nsh1 = cs1[c1].size();

        for (auto c2 = 0ul; c2 < csize2; ++c2) {
          auto sh3 = 0;
          const auto nsh2 = cs2[c2].size();

          for (auto c3 = 0ul; c3 < csize3; ++c3) {
            const auto nsh3 = cs3[c3].size();

            auto task_f = [=](float *out, float scaling_factor) {
              auto &Qab = Qab_->Q();
              auto &Qcd = Qcd_->Q();

              float norm = 0.0;
              for (auto a = sh0; a < nsh0 + sh0; ++a) {
                for (auto b = sh1; b < nsh1 + sh1; ++b) {
                  const auto ab = Qab(a, b);

                  if (ab * Qcd_->operator()() >= thresh2_) {
                    for (auto c = sh2; c < nsh2 + sh2; ++c) {
                      if (ab * Qcd_->max_in_row(c) >= thresh2_) {
                        for (auto d = sh3; d < nsh3 + sh3; ++d) {
                          const auto val = ab * Qcd(c, d);

                          norm += val;
                        }
                      }
                    }
                  }
                }
              }

              *out = scaling_factor * std::sqrt(norm);
            };

            if (plist.is_canonical(c0,c1,c2)) {
              const float multiplicity = static_cast<float>(plist.multiplicity(c0,c1,c2,c3));
              world.taskq.add(task_f, &norms(c0, c1, c2, c3), multiplicity);
            }

            sh3 += nsh3;
          }
          sh2 += nsh2;
        }
        sh1 += nsh1;
      }
      sh0 += nsh0;
    } // End estimate
  } else {
    norms = Screener::norm_estimate(world, bs_array, plist);
  }
  world.gop.fence();


  return norms;
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
