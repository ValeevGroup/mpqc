
#ifndef _chemistry_qc_wfn_wfn_h
#define _chemistry_qc_wfn_wfn_h

#include <stdio.h>
#include <util/misc/compute.h>
#include <math/scmat/matrix.h>
#include <math/topology/point.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/basis.h>

class Wavefunction: public MolecularEnergy
{
#   define CLASSNAME Wavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 private:
    RefSCDimension _basisdim;
    ResultRefSCMatrix _natural_orbitals;
    ResultRefDiagSCMatrix _natural_density;

    double* bs_values;
    double* bsg_values;
    RefGaussianBasisSet _gbs;
 public:
    Wavefunction(KeyVal&);
    Wavefunction(StateIn&);
    virtual ~Wavefunction();
    void save_data_state(StateOut&);

    void print(SCostream& =SCostream::cout);
    double density(cart_point&);
    double density_gradient(cart_point&,double*);
    double natural_orbital(cart_point& r, int iorb);
    double natural_orbital_density(cart_point& r, int orb, double* orbval = 0);
    double orbital(cart_point& r, int iorb, RefSCMatrix& orbs);
    double orbital_density(cart_point& r,
                           int iorb,
                           RefSCMatrix& orbs,
                           double* orbval = 0);

    virtual RefSymmSCMatrix density() = 0;
    virtual RefSCMatrix natural_orbitals();
    virtual RefDiagSCMatrix natural_density();
    RefGaussianBasisSet basis();
};

#endif
