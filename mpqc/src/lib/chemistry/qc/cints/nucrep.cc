
#include <chemistry/molecule/localdef.h>
#include <chemistry/qc/cints/cints.h>

double
cints_nuclear_repulsion_energy(const RefMolecule& mol_)
{
  int i, j;
  double r, e=0;

  Molecule& mol = *mol_.pointer();
  
  for (i=1; i < mol.natom(); i++) {
    AtomicCenter& ai = mol.atom(i);
    double Zi = ai.element().charge();
    
    for (j=0; j < i; j++) {
      AtomicCenter& aj = mol.atom(j);

      r = dist(ai.point(), aj.point());
      e += Zi * aj.element().charge() / dist(ai.point(), aj.point());
    }
  }

  return e;
}

double
cints_nuclear_repulsion_energy(const RefGaussianBasisSet& gbs)
{
  return cints_nuclear_repulsion_energy(gbs->molecule());
}

