
#include "mpqc/chemistry/molecule/energy.h"

namespace mpqc{

Energy::Energy(const KeyVal &kv) {
  if (kv.exists_class("atoms")) {
    atoms_ = kv.class_ptr<Molecule>("atoms");
  }
}

Energy::~Energy() = default;

}  // namespace mpqc
