

#include "mpqc/math/linalg/sqrt_inv.h"

namespace mpqc {
namespace math {

template<typename Array>
Array invert(Array const &S) {
  auto evg0 = tcc_time::now();
  auto spectral_range = pure::eval_guess(S);
  auto evg1 = tcc_time::now();
  auto eval_time = tcc_time::duration_in_s(evg0, evg1);
  if (S.world().rank() == 0) {
    std::cout << "\tEigenvalue estimation time = " << eval_time << " s\n";
  }

  const auto max_eval = spectral_range[1];
  const auto min_eval = std::max(0.0, spectral_range[0]);
  std::cout << "Min eval: " << min_eval << " Max eval: " << max_eval
            << std::endl;
  auto scale = 1.0 / (max_eval + min_eval);

  Array X, Inv;
  X("i,j") = scale * S("i,j");
  Inv("i,j") = X("i,k") * S("k,j");

  auto iter = 0;
  double trace_ideal = S.trange().tiles_range().extent()[0];
  double trace_real = 0.0;
  while (iter < 1000 && std::abs(trace_real - trace_ideal) >= 1e-10) {
    X("i,j") = 2 * X("i,j") - X("i,k") * S("k,l") * X("l,j");
    Inv("i,j") = X("i,k") * S("k,j");
    trace_real = Inv("i,j").trace();
    std::cout << "Error(" << iter
              << ") = " << std::abs(trace_real - trace_ideal) << std::endl;
    ++iter;
  }

  X.world().gop.fence();
  return X;
}

}  // namespace math
}  // namespace mpqc
