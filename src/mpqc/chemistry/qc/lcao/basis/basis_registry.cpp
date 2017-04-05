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

  ExEnv::out0() << "\nConstructing Basis Registry: \n";
  auto check = [](const std::shared_ptr<gaussian::AtomicBasis>& basis,
                  const char* basis_keyword) {
    if (!basis)
      throw InputError(
          (std::string("OrbitalBasisRegistry: keyword \"") + basis_keyword +
           "\" given, but not able to " + "construct AtomicBasis")
              .c_str(),
          __FILE__, __LINE__, basis_keyword);
  };

  auto basis = kv.class_ptr<gaussian::AtomicBasis>("basis");
  if (!basis)
    throw InputError("OrbitalBasisRegistry: missing basis", __FILE__, __LINE__,
                     "basis");
  this->add(OrbitalIndex(L"μ"), basis);
  ::mpqc::detail::parallel_print_range_info(
      world, basis->create_trange1(),
      "OBS Basis = " + kv.value<std::string>("basis:name"));

  if (kv.exists("df_basis")) {
    auto df_basis = kv.class_ptr<gaussian::AtomicBasis>("df_basis");
    check(df_basis, "df_basis");
    ::mpqc::detail::parallel_print_range_info(
        world, df_basis->create_trange1(),
        "DF Basis = " + kv.value<std::string>("df_basis:name"));
    this->add(OrbitalIndex(L"Κ"), df_basis);
  }

  if (kv.exists("aux_basis")) {
    auto aux_basis = kv.class_ptr<gaussian::AtomicBasis>("aux_basis");
    check(aux_basis, "aux_basis");
    this->add(OrbitalIndex(L"α"), aux_basis);
    ::mpqc::detail::parallel_print_range_info(
        world, aux_basis->create_trange1(),
        "AUX Basis = " + kv.value<std::string>("aux_basis:name"));
  }

  if (kv.exists("vir_basis")) {
    auto vir_basis = kv.class_ptr<gaussian::AtomicBasis>("vir_basis");
    check(vir_basis, "vir_basis");
    this->add(OrbitalIndex(L"Α"), vir_basis);
    ::mpqc::detail::parallel_print_range_info(
        world, vir_basis->create_trange1(),
        "Virtual Basis = " + kv.value<std::string>("vir_basis:name"));
  }
}

template <>
void gaussian::OrbitalBasisRegistry::clear() {
  // purge bases that do not update their state automatically, i.e.
  // non-AtomicBasis bases
  this->purge_if([](const value_type& val) {
    return !static_cast<bool>(
        std::dynamic_pointer_cast<gaussian::AtomicBasis>(val.second));
  });
}

}  // namespace lcao
}  // namespace mpqc
