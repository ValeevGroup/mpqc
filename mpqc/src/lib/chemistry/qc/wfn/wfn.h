
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
#include <chemistry/qc/basis/integral.h>

class Wavefunction: public MolecularEnergy {
#   define CLASSNAME Wavefunction
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefSCDimension basisdim_;
    RefSCMatrixKit basiskit_;
    
    ResultRefSymmSCMatrix overlap_;
    ResultRefSymmSCMatrix hcore_;
    ResultRefSCMatrix natural_orbitals_;
    ResultRefDiagSCMatrix natural_density_;

    double * bs_values;
    double * bsg_values;

    RefGaussianBasisSet gbs_;
    RefIntegral integral_;
    
  public:
    Wavefunction(StateIn&);
    Wavefunction(const RefKeyVal&);
    virtual ~Wavefunction();

    void save_data_state(StateOut&);

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
    virtual RefSymmSCMatrix core_hamiltonian();

    RefSCDimension basis_dimension();
    RefSCMatrixKit basis_matrixkit();
    RefGaussianBasisSet basis();
    RefIntegral integral();

    // returns a matrix which transforms AO's to orthogonal AO's
    // can be overridden, but defaults to S^-1/2
    virtual RefSymmSCMatrix ao_to_orthog_ao();

    void print(ostream& = cout);
};
SavableState_REF_dec(Wavefunction);

#endif
