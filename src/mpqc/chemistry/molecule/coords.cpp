
#include "mpqc/chemistry/molecule/coords.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("CartMolecularCoordinates", mpqc::CartMolecularCoordinates);

using namespace mpqc;

MolecularCoordinates::MolecularCoordinates(const std::shared_ptr<Molecule>& mol)
    : molecule_(mol) {}

MolecularCoordinates::MolecularCoordinates(const KeyVal& kv) {
  molecule_ = kv.class_ptr<Molecule>("molecule");
  if (!molecule_)
    throw InputError(
        "MolecularCoordinates KeyVal ctor did not receive a Molecule object",
        __FILE__, __LINE__, "molecule");
}

MolecularCoordinates::~MolecularCoordinates() {}

size_t MolecularCoordinates::nconstrained() const { return 0; }

CartMolecularCoordinates::CartMolecularCoordinates(
    const std::shared_ptr<Molecule>& mol)
    : MolecularCoordinates(mol) {}

CartMolecularCoordinates::CartMolecularCoordinates(
    const KeyVal& kv) : MolecularCoordinates(kv) {}

CartMolecularCoordinates::~CartMolecularCoordinates() {}

size_t CartMolecularCoordinates::size() const {
  return 3 * this->molecule()->atoms().size();
}
