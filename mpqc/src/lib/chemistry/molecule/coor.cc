
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

MolecularCoor::MolecularCoor(Molecule&mol):
  _mol(mol)
{
}
MolecularCoor::~MolecularCoor()
{
}

ProjectedCartesian::ProjectedCartesian(Molecule&mol):
  MolecularCoor(mol)
{
  int i;
  ColumnVector v[6];
  int nproj = 3;

  int natom = _mol.natom();
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

#if 0
NonRedundantIntCo::NonRedundantIntCo(Molecule&mol)
  : InternalCoordinate(mol), symm(0), nint(0)
{
  symm = Geom_form_symm(_mol);

  if(!symm)
    err_quit("NonRedundantIntCo(Molecule&mol): could not form symm");

  Geom_calc_simples(symm,_mol);

  for(SymmCoListIter p=symm; p; p++,nint++) ;
}

NonRedundantIntCo::~NonRedundantIntCo()
{
  if(symm) delete symm;
}

void NonRedundantIntCo::to_cartesian(ColumnVector&cartesian,
                                     ColumnVector&projected)
{
}

void NonRedundantIntCo::to_internal(ColumnVector&internal,
                                    ColumnVector&cartesian)
{
  DMatrix bmat = Geom_make_bmat(symm,_mol);

  DMatrix ginv = (bmat*bmat.transpose()).inverse();

  Dmatrix ctmp(bmat.ncol());
  for(int i=0; i < ctmp.dim(); i++) ctmp[i] = cartesian(i);
  
  DMatrix itmp = ginv*bmat*ctmp;

  for(i=0; i < itmp.dim(); i++) internal(i) = itmp[i];
}
#endif
