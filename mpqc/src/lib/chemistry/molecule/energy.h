
#ifndef _chemistry_molecule_energy_h
#define _chemistry_molecule_energy_h

#ifdef __GNUC__
#pragma interface
#endif

#include <iostream.h>

#include <math/optimize/function.h>
#include <math/optimize/conv.h>
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

    RefSCVector cartesian_gradient_;
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

    RefSCVector get_cartesian_x();
    RefSCVector get_cartesian_gradient();

    virtual void print(ostream& = cout);
};
SavableState_REF_dec(MolecularEnergy);

class MolEnergyConvergence: public Convergence {
#   define CLASSNAME MolEnergyConvergence
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int cartesian_;

    void set_defaults();
  public:
    // Standard constructors and destructor.
    MolEnergyConvergence();
    MolEnergyConvergence(StateIn&);
    MolEnergyConvergence(const RefKeyVal&);
    virtual ~MolEnergyConvergence();

    void save_data_state(StateOut&);

    // Set the current gradient and position information.  These
    //will possibly grab the cartesian infomation if the Function
    //is a MolecularEnergy.
    void get_grad(const RefFunction &);
    void get_x(const RefFunction &);
    void get_nextx(const RefFunction &);

    // Return nonzero if the optimization has converged.
    int converged();
};
SavableState_REF_dec(MolEnergyConvergence);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
