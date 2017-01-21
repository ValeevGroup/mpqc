//
// Created by Chong Peng on 9/9/16.
//

#include "basis_registry.h"
#include "mpqc/math/external/tiledarray/array_info.h"

namespace mpqc {
namespace lcao {

template <>
OrbitalRegistry<gaussian::Basis>::OrbitalRegistry(const KeyVal& kv)
    : Registry<OrbitalIndex, gaussian::Basis>() {
  auto& world = *kv.value<madness::World*>("$:world");

  auto basis = kv.class_ptr<gaussian::Basis>("basis");
  assert(basis != nullptr);
  this->add(OrbitalIndex(L"μ"), *basis);
  detail::parallel_print_range_info(world, basis->create_trange1(),
                                    "OBS Basis");

  if (kv.exists("df_basis")) {
    auto df_basis = kv.class_ptr<gaussian::Basis>("df_basis");
    assert(df_basis != nullptr);
    detail::parallel_print_range_info(world, df_basis->create_trange1(),
                                      "DF Basis");
    this->add(OrbitalIndex(L"Κ"), *df_basis);
  }

  if (kv.exists("aux_basis")) {
    auto aux_basis = kv.class_ptr<gaussian::Basis>("aux_basis");
    assert(aux_basis != nullptr);
    this->add(OrbitalIndex(L"α"), *aux_basis);
    detail::parallel_print_range_info(world, aux_basis->create_trange1(),
                                      "AUX Basis");
  }

  if (kv.exists("vir_basis")) {
    auto vir_basis = kv.class_ptr<gaussian::Basis>("vir_basis");
    assert(vir_basis != nullptr);
    this->add(OrbitalIndex(L"Α"), *vir_basis);
    detail::parallel_print_range_info(world, vir_basis->create_trange1(),
                                      "Virtual Basis");
  }
}

}  // namespace lcao
}  // namespace mpqc
