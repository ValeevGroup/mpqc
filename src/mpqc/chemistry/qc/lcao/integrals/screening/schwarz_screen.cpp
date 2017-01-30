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

bool SchwarzScreen::validate(int64_t a, int64_t b, int64_t c, 
                             double const* buf, int64_t size) const {
  // First compute relevant values
  auto norm = 0.0;
  for (auto i = 0; i < size; ++i) { 
    const auto val = buf[i]; // May not hold due to diagonal estimators
    norm = std::max(norm, std::abs(val));
  }
  const auto norm2 = norm * norm;

  const auto est0 = estimate(a);
  const auto est1 = estimate(a, b);
  const auto est2 = estimate(a, b, c);

  auto failed = false;
  if (est0 < norm2 || est1 < norm2 || est2 < norm2) {
    failed = true;
  }
  if (failed && thresh2_ < norm2) {
    ExEnv::out0() << "\nSize = " << size << std::endl;
    ExEnv::out0() << "Qab = " << Qab()(a) << ", Qcd = " << Qcd()(b, c)
              << std::endl;
    ExEnv::out0() << "Integral screening for " << a << ", " << b << ", " << c
              << "\n\twith ests " << est0 << ", " << est1 << ", "
              << est2 << "\n\twas inaccurate for norm " << norm
              << std::endl;
    ExEnv::out0() << "Integral values:\n\t";
    for (auto i = 0; i < size; ++i) {
      ExEnv::out0() << buf[i] << " ";
    }
    ExEnv::out0() << std::endl;
    throw;
  }

  return failed;
}

bool SchwarzScreen::validate(int64_t a, int64_t b, int64_t c, int64_t d,
                             double const* buf, int64_t size) const {
  // First compute relevant values
  auto norm = 0.0;
  for (auto i = 0; i < size; ++i) { // compute the inf norm of the shell quartet
    const auto val = buf[i];
    norm = std::max(norm, std::abs(val));
  }
  const auto norm2 = norm * norm; // We compute squared estimates

  // Grab all the used estimates
  const auto est0 = estimate(a);
  const auto est1 = estimate(a, b);
  const auto est2 = estimate(a, b, c);
  const auto est3 = estimate(a, b, c, d);

  auto failed = false; // Assume good for now
  if (est0 < norm2 || est1 < norm2 || est2 < norm2 || est3 < norm2) {
    failed = true; // If not good then failed
  }
  if (failed && thresh2_ < norm2) {
    ExEnv::out0() << "\nSize = " << size << std::endl;
    ExEnv::out0() << "Qab = " << Qab()(a, b) << ", Qcd = " << Qcd()(c, d)
              << std::endl;
    ExEnv::out0() << "Integral screening for " << a << ", " << b << ", " << c
              << ", " << d << "\n\twith ests " << est0 << ", " << est1 << ", "
              << est2 << ", " << est3 << "\n\twas inaccurate for norm " << norm
              << std::endl;
    ExEnv::out0() << "Integral values:\n\t";
    for (auto i = 0; i < size; ++i) {
      ExEnv::out0() << buf[i] << " ";
    }
    std::cout << std::endl;
    throw; // TODO figure out what this should do
  }

  return failed;
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
