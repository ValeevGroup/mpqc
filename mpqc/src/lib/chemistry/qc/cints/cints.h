
#ifndef _chemistry_qc_cints_h
#define _chemistry_qc_cints_h

#include <chemistry/qc/basis/gaussbas.h>

#define MAXAM 5
#define MAXCLASS 1000

double cints_nuclear_repulsion_energy(const RefMolecule& mol_);
double cints_nuclear_repulsion_energy(const RefGaussianBasisSet& gbs);
void cints_shell_overlap(const RefGaussianBasisSet& bs1,
                         const RefGaussianBasisSet& bs2,
                         double *buf, int i, int j);

#endif
