
#ifndef _chemistry_molecule_taylor_h
#define _chemistry_molecule_taylor_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/molecule/energy.h>
#include <chemistry/molecule/coor.h>

// the molecular energy as a taylor expansion
class TaylorMolecularEnergy: public MolecularEnergy {
#   define CLASSNAME TaylorMolecularEnergy
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    void compute_energy(double&);
    
    // the coordinates
    RefSetIntCoor coordinates_;

    // The force constants (only the unique ones are given) to arbitrary
    // order.  If nonunique force constants are put here, then the answer
    // will be wrong
    ArrayArrayint force_constant_index_;
    Arraydouble force_constant_value_;
    
    // the dimension of coordinates_;
    RefSCDimension dim_;

    // the expansion point
    RefSCVector expansion_point_;

    // the energy at the expansion point
    double e0_;
  public:
    TaylorMolecularEnergy(KeyVal&);
    TaylorMolecularEnergy(StateIn&);
    ~TaylorMolecularEnergy();
    void save_data_state(StateOut&);
    void print(SCostream& = SCostream::cout);
    void compute();
};

#endif
