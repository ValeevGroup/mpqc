
#ifndef _chemistry_qc_wfn_obwfn_h
#define _chemistry_qc_wfn_obwfn_h

#include <chemistry/qc/wfn/wfn.h>

class OneBodyWavefunction: public Wavefunction
{
 private:
  void init(KeyVal&keyval);
  int _have_density;
  SymmetricMatrix _density;
 public:
  OneBodyWavefunction(KeyVal&,Molecule&,GaussianBasisSet&);
  OneBodyWavefunction(KeyVal&,Molecule&,GaussianBasisSet&,MolecularCoor&);
  ~OneBodyWavefunction();
  void x_changed();

  virtual const Matrix& eigenvectors() = 0;
  virtual double occupation(int vectornum) = 0;
  double orbital(cart_point& r, int iorb);
  double orbital_density(cart_point& r, int iorb, double* orbval = 0);

  const SymmetricMatrix& density();
};

#endif
