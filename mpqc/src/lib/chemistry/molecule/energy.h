
#ifndef _chemistry_molecule_energy_h
#define _chemistry_molecule_energy_h

#ifdef __GNUC__
#pragma interface
#endif

extern "C" {
#include <stdio.h>
}

#include <math/optimize/nlp.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>

class MolecularEnergy: public NLP2 {
#   define CLASSNAME MolecularEnergy
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  private:
    RefMolecularCoor _mc;

  protected:
    AccResultdouble& _energy;

    RefSCDimension _moldim; // the number of cartesian variables
    RefMolecule _mol;

    void failure(const char *);

    virtual void set_energy(double);
    // These are passed gradients and hessian in cartesian coordinates.
    // The _gradient and _hessian in internal coordinates are computed.
    virtual void set_gradient(RefSCVector&);
    virtual void set_hessian(RefSymmSCMatrix&);

    void x_to_molecule();
    void molecule_to_x();
  public:
    MolecularEnergy(const RefKeyVal&);
    MolecularEnergy(StateIn&);
    MolecularEnergy(RefMolecule&);
    MolecularEnergy(RefMolecule&,RefMolecularCoor&);
    ~MolecularEnergy();
    void save_data_state(StateOut&);

    virtual double energy();
    virtual RefMolecule molecule();
    void guess_hessian(RefSymmSCMatrix&);
    RefSymmSCMatrix inverse_hessian(RefSymmSCMatrix&);

    void set_x(const RefSCVector&);

    virtual void print(SCostream& =SCostream::cout);
};
SavableState_REF_dec(MolecularEnergy);

#endif
