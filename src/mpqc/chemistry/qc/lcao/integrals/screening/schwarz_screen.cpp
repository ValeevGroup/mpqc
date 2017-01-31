#include "mpqc/chemistry/qc/lcao/integrals/screening/schwarz_screen.h"
#include <mpqc/util/misc/exenv.h>

namespace mpqc {
namespace lcao {
namespace gaussian {

int64_t Qmatrix::f2s(Qmatrix::f2s_map const& map, int64_t ind) const {
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

boost::optional<double> SchwarzScreen::estimate(int64_t a,
                                                       int64_t b) const {
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

boost::optional<double> SchwarzScreen::estimate(int64_t a, int64_t b,
                                                       int64_t c,
                                                       int64_t d) const {
  assert(!Qab_->is_aux_Q());
  return Qab()(a, b) * Qcd()(c, d);
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
