
#ifndef _chemistry_molecule_energy_h
#define _chemistry_molecule_energy_h

#ifdef __GNUC__
#pragma interface
#endif

extern "C" {
#include <stdio.h>
}

#include <math/optimize/function.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

class MolecularEnergy: public Function {
#   define CLASSNAME MolecularEnergy
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefSCDimension moldim_; // the number of cartesian variables
    RefMolecularCoor mc_;
    RefMolecule mol_;

  protected:
    void failure(const char *);

    // this is just a wrapper around set_value()
    virtual void set_energy(double);

    // These are passed gradients and hessian in cartesian coordinates.
    // The _gradient and _hessian in internal coordinates are computed.
    virtual void set_gradient(RefSCVector&);
    virtual void set_hessian(RefSymmSCMatrix&);

    void x_to_molecule();
    void molecule_to_x();

  public:
    MolecularEnergy(const MolecularEnergy&);
    MolecularEnergy(const RefKeyVal&);
    MolecularEnergy(StateIn&);
    ~MolecularEnergy();

    void save_data_state(StateOut&);

    MolecularEnergy & operator=(const MolecularEnergy&);
    
    // a wrapper around value();
    virtual double energy();

    virtual RefMolecule molecule();
    virtual RefSCDimension moldim();
    
    void guess_hessian(RefSymmSCMatrix&);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    void set_x(const RefSCVector&);

    virtual void print(ostream& = cout);
};
SavableState_REF_dec(MolecularEnergy);

#endif
