
#ifndef _chemistry_qc_mpqc_h
#define _chemistry_qc_mpqc_h

#include <stdio.h>
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/molecule/energy.h>
#include <math/newmat7/newmat.h>
#include <math/topology/point.h>
#include <chemistry/qc/basis/basis.h>
class Molecule;
class MolecularCoor;
class ColumnVector;
class SymmetricMatrix;
class Matrix;

class MPQC: public OneBodyWavefunction
{
 private:
  int _have_exchange_energy;
  int _have_eigenvectors;
  int _have_eigenvectors_on_disk;
  int _have_old_eigenvectors_on_disk;
  int _do_exchange_energy;
  int _do_eigenvectors;
  int _do_eigenvectors_on_disk;
  int _do_reuse_old_eigenvectors;
  double* bs_values;
  double* bsg_values;
  int _nocc;
  int _maxiter;
  double _exchange_energy;
  void init(KeyVal&);
 protected:
  Matrix _eigenvectors;
  void compute();
  void read_vector(char *fname,int n_basis, Matrix &scf_vector);
 public:
  MPQC(KeyVal&,Molecule&,GaussianBasisSet&);
  MPQC(KeyVal&,Molecule&,GaussianBasisSet&,MolecularCoor&);
  virtual ~MPQC();

  void x_changed();

  void print(FILE* =stdout);
  int do_exchange_energy(int);
  int do_eigenvectors(int);
  int do_eigenvectors_on_disk(int);
  int do_reuse_old_eigenvectors(int);
  const Matrix& eigenvectors();
  double exchange_energy();
  void maxiter(int m);

  double occupation(int vectornum);
};

#endif
