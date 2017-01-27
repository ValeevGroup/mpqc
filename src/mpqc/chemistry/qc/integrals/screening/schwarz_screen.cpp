#include "mpqc/chemistry/qc/integrals/screening/schwarz_screen.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

/* Qmatrix::Qmatrix(RowMatrixXd Q, std::unordered_map<int64_t, int64_t> map) */
/*     : Q_(std::move(Q)), */
/*       func_to_shell_map_(std::move(map)), */
/*       max_elem_in_row_(Eigen::VectorXd(Q_.rows())), */
/*       max_elem_in_Q_(0.0) */

/* { */
/*   const auto nrows = Q_.rows(); */
/*   for (auto i = 0; i < nrows; ++i) { */
/*     const auto max_elem = Q_.row(i).cwiseAbs().maxCoeff(); */
/*     max_elem_in_Q_ = std::max(max_elem, max_elem_in_Q_); */
/*     max_elem_in_row_(i) = max_elem; */
/*   } */
/* } */

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

bool Qmatrix::is_aux_Q() const {
  return Q_.rows() == 1;
}

SchwarzScreen::SchwarzScreen(std::shared_ptr<Qmatrix> Qab,
                               std::shared_ptr<Qmatrix> Qcd, double thresh)
    : Screener(), Qab_(std::move(Qab)), Qcd_(std::move(Qcd)), thresh_(thresh) {
  if (Qab_->is_aux_Q()) {
    Q_ab_is_aux_ = true;
  }
}

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

inline boost::optional<double> SchwarzScreen::estimate(int64_t a) const {
  return Qab()(a) * Qcd()();
}

inline boost::optional<double> SchwarzScreen::estimate(int64_t a,
                                                        int64_t b) const {
  if (Q_ab_is_aux_) {
    return Qab()(a) * Qcd()(b);
  } else {
    return Qab()(a, b) * Qcd()();
  }
}

inline boost::optional<double> SchwarzScreen::estimate(int64_t a, int64_t b,
                                                        int64_t c) const {
  if(Q_ab_is_aux_){
    return Qab()(a) * Qcd()(b, c);
  } else {
    return Qab()(a,b) * Qcd()(c);
  }
}

inline boost::optional<double> SchwarzScreen::estimate(int64_t a, int64_t b,
                                                        int64_t c,
                                                        int64_t d) const {
  assert(!Q_ab_is_aux_);
  return Qab()(a, b) * Qcd()(c, d);
}

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
