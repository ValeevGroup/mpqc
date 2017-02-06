//
// Created by Chong Peng on 9/9/16.
//

#include "basis_registry.h"
#include "mpqc/math/external/tiledarray/array_info.h"

namespace mpqc {
namespace lcao {

template <>
gaussian::OrbitalBasisRegistry::OrbitalRegistry(const KeyVal& kv)
    : Registry<OrbitalIndex, std::shared_ptr<gaussian::Basis>>() {
  auto& world = *kv.value<madness::World*>("$:world");

  auto basis = kv.class_ptr<gaussian::AtomicBasis>("basis");
  assert(basis != nullptr);
  this->add(OrbitalIndex(L"μ"), basis);
  ::mpqc::detail::parallel_print_range_info(world, basis->create_trange1(),
                                            "OBS Basis");

  if (kv.exists("df_basis")) {
    auto df_basis = kv.class_ptr<gaussian::AtomicBasis>("df_basis");
    assert(df_basis != nullptr);
    ::mpqc::detail::parallel_print_range_info(
                         world, df_basis->create_trange1(), "DF Basis");
    this->add(OrbitalIndex(L"Κ"), df_basis);
  }

  if (kv.exists("aux_basis")) {
    auto aux_basis = kv.class_ptr<gaussian::AtomicBasis>("aux_basis");
    assert(aux_basis != nullptr);
    this->add(OrbitalIndex(L"α"), aux_basis);
    ::mpqc::detail::parallel_print_range_info(
                         world, aux_basis->create_trange1(), "AUX Basis");
  }

  if (kv.exists("vir_basis")) {
    auto vir_basis = kv.class_ptr<gaussian::AtomicBasis>("vir_basis");
    assert(vir_basis != nullptr);
    this->add(OrbitalIndex(L"Α"), vir_basis);
    ::mpqc::detail::parallel_print_range_info(
                         world, vir_basis->create_trange1(), "Virtual Basis");
  }
}

template <>
void gaussian::OrbitalBasisRegistry::clear() {
  // purge bases that do not update their state automatically, i.e. non-AtomicBasis bases
  this->purge_if([](const value_type& val) -> bool {
    return !static_cast<bool>(
        std::dynamic_pointer_cast<gaussian::AtomicBasis>(val.second));
  });
}

}  // namespace lcao
}  // namespace mpqc
