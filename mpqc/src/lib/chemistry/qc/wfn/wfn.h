
#ifndef _chemistry_qc_wfn_wfn_h
#define _chemistry_qc_wfn_wfn_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <util/misc/compute.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/energy.h>
#include <chemistry/qc/basis/basis.h>

class Wavefunction: public MolecularEnergy
{
#   define CLASSNAME Wavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
 private:
    RefSCDimension _basisdim;
    ResultRefSymmSCMatrix _overlap;
    ResultRefSCMatrix _natural_orbitals;
    ResultRefDiagSCMatrix _natural_density;

    double* bs_values;
    double* bsg_values;
    RefGaussianBasisSet _gbs;
 public:
    Wavefunction(const Wavefunction&);
    Wavefunction(const RefKeyVal&);
    Wavefunction(StateIn&);
    virtual ~Wavefunction();

    Wavefunction & operator=(const Wavefunction&);
    
    void save_data_state(StateOut&);

    void print(ostream& = cout);
    double density(const SCVector3&);
    double density_gradient(const SCVector3&,double*);
    double natural_orbital(const SCVector3& r, int iorb);
    double natural_orbital_density(const SCVector3& r,
                                   int orb, double* orbval = 0);
    double orbital(const SCVector3& r, int iorb, const RefSCMatrix& orbs);
    double orbital_density(const SCVector3& r,
                           int iorb,
                           const RefSCMatrix& orbs,
                           double* orbval = 0);

    virtual RefSymmSCMatrix density() = 0;
    virtual RefSCMatrix natural_orbitals();
    virtual RefDiagSCMatrix natural_density();
    virtual RefSymmSCMatrix overlap();
    RefGaussianBasisSet basis();

    RefSCDimension basis_dimension() { return _basisdim; }
};
SavableState_REF_dec(Wavefunction);

#endif
