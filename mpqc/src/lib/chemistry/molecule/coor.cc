
extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
}

#include <util/unix/cct_cprot.h>

#include <math/newmat7/newmat.h>
#include <math/nihmatrix/nihmatrix.h>
#include "molecule.h"
#include "coor.h"
#include "simple.h"
#include "symm.h"
#include "symmQCList.h"

MolecularCoor::MolecularCoor(RefMolecule mol):
  _mol(mol)
{
}
MolecularCoor::~MolecularCoor()
{
}

ProjectedCartesian::ProjectedCartesian(RefMolecule&mol):
  MolecularCoor(mol)
{
  int i;
  ColumnVector v[6];
  int nproj = 3;

  int natom = _mol->natom();
  ncart = 3*natom;
  nint = ncart - nproj;
  for (i=0; i<nproj; i++) v[i].ReDimension(ncart);
  double one_over_srnatom = 1.0/sqrt(natom);
  for (i=0; i<natom; i++) {
      v[0](i*3+1) = one_over_srnatom;
      v[0](i*3+2) = 0.0;
      v[0](i*3+3) = 0.0;
      v[1](i*3+1) = 0.0;
      v[1](i*3+2) = one_over_srnatom;
      v[1](i*3+3) = 0.0;
      v[2](i*3+1) = 0.0;
      v[2](i*3+2) = 0.0;
      v[2](i*3+3) = one_over_srnatom;
    }

  x.ReDimension(ncart,nint);
  // initial x to be a partial unit matrix
 for (i=1; i<=ncart; i++) {
     for (int j=1; j<=nint; j++) {
	  if (j==i) x(i,j) = 1.0;
	  else x(i,j) = 0.0;
	}
   }

  Print(v[0]);
  Print(v[1]);
  Print(v[2]);
  Print(x);

  // project out the v vectors from each of the rows of x
 // the v must be normalized
 for (i=1; i<nint; i++) {
     for (int j=0; j<nproj; j++) {
	 int k;
	 double dot = 0.0;
	 for (k=1; k<=ncart; k++) dot += x(k,i)*v[j](k);
	 for (k=1; k<=ncart; k++) {
	     x(k,i) = x(k,i) - dot*v[j](k);
	   }
       }
   }

  Print(x);

}

ProjectedCartesian::~ProjectedCartesian()
{
}

int ProjectedCartesian::dim()
{
  return nint;
}

void ProjectedCartesian::to_cartesian(ColumnVector&cartesian,
				      ColumnVector&projected)
{
  for (int i=1; i<=ncart; i++) {
      cartesian(i) = 0.0;
      for (int j=1; j<=nint; j++) {
	  cartesian(i) += projected(j)*x(i,j);
	}
    }
}

void ProjectedCartesian::to_internal(ColumnVector&projected,
				     ColumnVector&cartesian)
{
  for (int i=1; i<=nint; i++) {
      projected(i) = 0.0;
      for (int j=1; j<=ncart; j++) {
	  projected(i) += x(j,i)*cartesian(j);
	}
    }
}

////////////////////////////////////////////////////////////////////////////////
