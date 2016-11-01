//
// Created by Chong Peng on 9/9/16.
//

#include "basis_registry.h"
#include "../../../../../utility/cc_utility.h"

namespace mpqc {

template <>
mpqc::OrbitalRegistry<mpqc::basis::Basis>::OrbitalRegistry(const KeyVal& kv)
    : Registry<OrbitalIndex, mpqc::basis::Basis>() {
  auto& world = *kv.value<madness::World*>("$:world");

  auto basis = kv.keyval("basis").class_ptr<basis::Basis>();
  assert(basis != nullptr);
  this->add(OrbitalIndex(L"μ"), *basis);
  cc::parallel_print_range_info(world, basis->create_trange1(), "OBS Basis");

  if (kv.exists("df_basis")) {
    auto df_basis = kv.keyval("df_basis").class_ptr<basis::Basis>();
    assert(df_basis != nullptr);
    cc::parallel_print_range_info(world, df_basis->create_trange1(),
                                  "DF Basis");
    this->add(OrbitalIndex(L"Κ"), *df_basis);
  }

  if (kv.exists("aux_basis")) {
    auto aux_basis = kv.keyval("aux_basis").class_ptr<basis::Basis>();
    assert(aux_basis != nullptr);
    this->add(OrbitalIndex(L"α"), *aux_basis);
    cc::parallel_print_range_info(world, aux_basis->create_trange1(),
                                  "AUX Basis");
  }

  if (kv.exists("vir_basis")) {
    auto vir_basis = kv.keyval("vir_basis").class_ptr<basis::Basis>();
    assert(vir_basis != nullptr);
    this->add(OrbitalIndex(L"Α"), *vir_basis);
    cc::parallel_print_range_info(world, vir_basis->create_trange1(),
                                  "Virtual Basis");
  }
}

}  // endof namespace mpqc
