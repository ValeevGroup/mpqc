
#include <stdio.h>
#include <iostream.h>
#include <math/newmat7/newmat.h>
#include <util/keyval/keyval.h>
#include <math/topology/point.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/gaussbas.h>
#include "wfn.h"

// Function for returning an orbital value at a point
double Wavefunction::natural_orbital(cart_point& r, int iorb)
{
  return orbital(r,iorb,natural_orbitals());
}

// Function for returning an orbital value at a point
double Wavefunction::natural_orbital_density(cart_point& r,
                                             int iorb,
                                             double* orbval)
{
  return orbital_density(r,iorb,natural_orbitals(),orbval);
}

// Function for returning an orbital value at a point
double Wavefunction::orbital(cart_point& r,
                             int iorb,
                             const Matrix& orbs)
{
    double *result;
    int nbasis = basis().nbasis();
    if (!bs_values) bs_values=new double[nbasis];

    // compute the basis set values
    basis().values(r,bs_values);
    
    // loop over basis functions
    double orb_value = 0;
    for (int i=0; i<nbasis; i++)
        orb_value += orbs.element(i,iorb)*bs_values[i];

    return orb_value;
}     

double Wavefunction::orbital_density(cart_point& r,
                                     int iorb,
                                     const Matrix& orbs,
                                     double* orbvalue)
{
  double tmp = orbital(r,iorb,orbs);
  if (orbvalue) *orbvalue = tmp;
  return 2.0 * tmp * tmp;
}
