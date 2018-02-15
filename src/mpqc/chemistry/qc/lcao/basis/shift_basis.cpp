#include "mpqc/chemistry/qc/lcao/basis/shift_basis.h"

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

std::shared_ptr<Basis> shift_basis_origin(const Basis &basis,
                                          const Vector3d &shift) {
  std::vector<ShellVec> vec_of_shells;
  for (auto shell_vec : basis.cluster_shells()) {
    ShellVec shells;
    for (auto shell : shell_vec) {
      std::array<double, 3> new_origin = {{shell.O[0] + shift(0),
                                           shell.O[1] + shift(1),
                                           shell.O[2] + shift(2)}};
      shell.move(new_origin);
      shells.push_back(shell);
    }
    vec_of_shells.push_back(shells);
  }
  Basis result(vec_of_shells);
  auto result_ptr = std::make_shared<Basis>(result);
  return result_ptr;
}

std::shared_ptr<Basis> shift_basis_origin(const Basis &basis,
                                          const Vector3d &shift_base,
                                          const Vector3i &nshift,
                                          const Vector3d &dcell,
                                          const bool is_half_range) {
  std::vector<ShellVec> vec_of_shells;

  using ::mpqc::detail::direct_ord_idx;
  using ::mpqc::detail::direct_vector;
  int64_t shift_size = 1 + direct_ord_idx(nshift, nshift);
  assert(shift_size > 0 && shift_size % 2 == 1);
  int64_t start = is_half_range ? ((shift_size - 1) / 2) : 0;

  for (auto idx_shift = start; idx_shift < shift_size; ++idx_shift) {
    Vector3d shift = direct_vector(idx_shift, nshift, dcell) + shift_base;

    for (auto shell_vec : basis.cluster_shells()) {
      ShellVec shells;
      for (auto shell : shell_vec) {
        std::array<double, 3> new_origin = {{shell.O[0] + shift(0),
                                             shell.O[1] + shift(1),
                                             shell.O[2] + shift(2)}};
        shell.move(new_origin);
        shells.push_back(shell);
      }
      vec_of_shells.push_back(shells);
    }
  }

  Basis result(vec_of_shells);
  auto result_ptr = std::make_shared<Basis>(result);
  return result_ptr;
}

}  // namespace detail
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
