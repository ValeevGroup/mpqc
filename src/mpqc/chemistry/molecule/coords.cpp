
#include "mpqc/chemistry/molecule/coords.h"
#include "mpqc/util/core/exception.h"
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

std::ostream& mpqc::operator<<(std::ostream& os, const MolecularCoordinates& coord) {
  coord.print(os);
  return os;
}

CartMolecularCoordinates::CartMolecularCoordinates(
    const std::shared_ptr<Molecule>& mol)
    : MolecularCoordinates(mol) {}

CartMolecularCoordinates::CartMolecularCoordinates(
    const KeyVal& kv) : MolecularCoordinates(kv) {}

CartMolecularCoordinates::~CartMolecularCoordinates() {}

size_t CartMolecularCoordinates::size() const {
  return 3 * this->molecule()->atoms().size();
}

std::shared_ptr<MolecularCoordinates> CartMolecularCoordinates::clone() const {
  return std::make_shared<CartMolecularCoordinates>(*this);
}

void CartMolecularCoordinates::displace(size_t ncoords, size_t* coords,
                                        double* displacements) {
  auto atoms = molecule()->atoms();
  for(size_t c=0; c!=ncoords; ++c) {
    const auto atom = coords[c] / 3;
    const auto xyz = coords[c] % 3;
    auto r = atoms[atom].center();
    r(xyz) += displacements[c];
    auto mass = atoms[atom].mass();
    auto Z = atoms[atom].atomic_number();
    auto charge = atoms[atom].charge();

    atoms[atom] = Atom(r, mass, Z, charge);
  }

  this->update_molecule(atoms);
}

void CartMolecularCoordinates::print(std::ostream& os) const {
  auto atoms = molecule()->atoms();
  os << indent << "CartMolecularCoordinates:" << std::endl << incindent;
  os << *molecule();
  os << decindent;
}
