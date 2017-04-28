#include "mpqc/chemistry/qc/lcao/scf/decomposed_rij.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

namespace detail {

std::shared_ptr<Basis> create_unit_basis() {
  std::vector<ShellVec> vec_of_shells;
  vec_of_shells.push_back({{Shell::unit()}});

  Basis result(vec_of_shells);

  return std::make_shared<Basis>(result);
}

}  // namespace detail

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
