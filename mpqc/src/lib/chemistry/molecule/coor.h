
#ifndef _chemistry_molecule_molcoor_h
#define _chemistry_molecule_molcoor_h

#include <stdio.h>
#include <math/newmat7/newmat.h>
class Molecule;

class MolecularCoor
{
 protected:
  RefMolecule _mol;
 public:
  MolecularCoor(RefMolecule mol);
  virtual ~MolecularCoor();
  virtual int dim() = 0;
  virtual void to_cartesian(ColumnVector&cartesian,ColumnVector&internal) = 0;
  virtual void to_internal(ColumnVector&internal,ColumnVector&cartesian) = 0;
};

class ProjectedCartesian: public MolecularCoor
{
 private:
  int nproj; // for linear 5, nonlinear 6
  Matrix x;
  int ncart;
  int nint;
 public:
  ProjectedCartesian(RefMolecule&mol);
  ~ProjectedCartesian();
  int dim();
  void to_cartesian(ColumnVector&cartesian,ColumnVector&projected);
  void to_internal(ColumnVector&projected,ColumnVector&cartesian);
};

/////////////////////////////////////////////////////////////////////

#if 0
class InternalMolecularCoor : public MolecularCoor
{
  protected:
    InternalCoorSet variable;
    InternalCoorSet fixed;
};
#endif

#endif
