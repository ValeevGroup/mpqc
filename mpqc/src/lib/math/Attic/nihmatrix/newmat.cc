
#include "nihmatrix.h"
#include <math/newmat7/newmat.h>


void Convert(DVector&in,ColumnVector&out)
{
  int dim = in.dim();
  out.ReDimension(dim);
  for (int i=0; i<dim; i++) out.element(i) = in[i];
}

void Convert(DMatrix&in,SymmetricMatrix&out)
{
  int dim = in.nrow();
  out.ReDimension(dim);
  for (int i=0; i<dim; i++) {
      for (int j=0; j<=i; j++) {
          out.element(i,j) = (in(i,j)+in(j,i))*0.5;
        }
    }
}

void Convert(ColumnVector&in,DVector&out)
{
  int dim = in.Nrows();
  out.resize(dim);
  for (int i=0; i<dim; i++) out[i] = in.element(i);
}

void Convert(SymmetricMatrix&in,DMatrix&out)
{
  int dim = in.Nrows();
  out.resize(dim,dim);
  for (int i=0; i<dim; i++) {
      for (int j=0; j<=i; j++) {
          out(i,j) = in.element(i,j);
        }
    }
}
