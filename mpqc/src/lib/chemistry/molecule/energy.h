
#ifndef _chemistry_molecule_energy_h
#define _chemistry_molecule_energy_h

extern "C" {
#include <stdio.h>
}

#include <math/opt/nlp.h>
class Molecule;
class ColumnVector;
class SymmetricMatrix;
class MolecularCoor;

class MolecularEnergy: public NLP2 {
  private:
    MolecularCoor* _mc;

    int _do_energy;
    int _do_gradient;
    int _do_hessian;

    int _have_energy;
    int _have_gradient;
    int _have_hessian;

    double& _energy;
    ColumnVector& _gradient;
    SymmetricMatrix& _hessian;
  protected:
    Molecule& _mol;

    void failure(const char *);

    virtual void compute() = 0;

    virtual void set_energy(double);
    // These are passed gradients and hessian in cartesian coordinates.
    // The _gradient and _hessian in internal coordinates are computed.
    virtual void set_gradient(ColumnVector&);
    virtual void set_hessian(SymmetricMatrix&);
  public:
    MolecularEnergy(Molecule&);
    MolecularEnergy(Molecule&,MolecularCoor&);
    ~MolecularEnergy();
    void X_to_molecule();
    void molecule_to_X();

    void Eval();
    double EvalF();
    ColumnVector EvalG();
    SymmetricMatrix EvalH();

    virtual void x_changed();

    virtual double energy();
    virtual const ColumnVector& gradient();
    virtual const SymmetricMatrix& hessian();

    int do_energy(int);
    int do_gradient(int);
    int do_hessian(int);
    int do_energy();
    int do_gradient();
    int do_hessian();

    inline Molecule& molecule() { return _mol; }
};

#endif
