//$$ newmat1.cxx   Matrix type functions

// Copyright (C) 1991,2,3: R B Davies


#include "include.h"

#include "newmat.h"



/************************* MatrixType functions *****************************/


MatrixType MatrixType::operator*(const MatrixType& mt) const
{
   int a = attribute & mt.attribute;
   if ((a & (Upper+Lower)) == (Upper+Lower)) return Dg;
   else return a & ~Symmetric;
}

MatrixType MatrixType::t() const
{
   int a = attribute & ~(Upper+Lower);
   if (attribute & Upper) a |= Lower;
   if (attribute & Lower) a |= Upper;
   return a;
}

MatrixType MatrixType::MultRHS() const
{
   if ((attribute & (Upper+Lower)) == (Upper+Lower)) return Dg;
   else return attribute & ~Symmetric;
}



MatrixType::operator char*() const
{
// make a string with the name of matrix with the given attributes
   switch (attribute)
   {
   case Valid:                              return "Rect ";
   case Valid+Symmetric:                    return "Sym  ";
   case Valid+Band:                         return "Band ";
   case Valid+Symmetric+Band:               return "SmBnd";
   case Valid+Upper:                        return "UT   ";
   case Valid+Upper+Lower:
   case Valid+Band+Upper+Lower:
   case Valid+Symmetric+Upper:
   case Valid+Symmetric+Band+Upper:
   case Valid+Symmetric+Lower:
   case Valid+Symmetric+Band+Lower:
   case Valid+Symmetric+Upper+Lower:
   case Valid+Symmetric+Band+Upper+Lower:   return "Diag ";
   case Valid+Band+Upper:                   return "UpBnd";
   case Valid+Lower:                        return "LT   ";
   case Valid+Band+Lower:                   return "LwBnd";
   default:
      if (!(attribute & Valid))             return "UnSp ";
      if (attribute & LUDeco)
         return (attribute & Band) ?     "BndLU" : "Crout";
                                            return "?????";
   }
}


GeneralMatrix* MatrixType::New(int nr, int nc, BaseMatrix* bm) const
{
// make a new matrix with the given attributes

   Tracer tr("New"); GeneralMatrix* gm;
   switch (attribute)
   {
   case Valid:
      if (nc==1) { gm = new ColumnVector(nr); break; }
      if (nr==1) { gm = new RowVector(nc); break; }
      gm = new Matrix(nr, nc); break;

   case Valid+Symmetric:
      gm = new SymmetricMatrix(nr); break;

   case Valid+Band:
      {
         MatrixBandWidth bw = bm->BandWidth();
         gm = new BandMatrix(nr,bw.lower,bw.upper); break;
      }

   case Valid+Symmetric+Band:
      gm = new SymmetricBandMatrix(nr,bm->BandWidth().lower); break;

   case Valid+Upper:
      gm = new UpperTriangularMatrix(nr); break;

   case Valid+Upper+Lower:
   case Valid+Band+Upper+Lower:
   case Valid+Symmetric+Upper:
   case Valid+Symmetric+Band+Upper:
   case Valid+Symmetric+Lower:
   case Valid+Symmetric+Band+Lower:
   case Valid+Symmetric+Upper+Lower:
   case Valid+Symmetric+Band+Upper+Lower:
      gm = new DiagonalMatrix(nr); break;

   case Valid+Band+Upper:
      gm = new UpperBandMatrix(nr,bm->BandWidth().upper); break;

   case Valid+Lower:
      gm = new LowerTriangularMatrix(nr); break;

   case Valid+Band+Lower:
      gm = new LowerBandMatrix(nr,bm->BandWidth().lower); break;

   default:
      Throw(ProgramException("Invalid matrix type"));
   }
   
   MatrixErrorNoSpace(gm); gm->Protect(); return gm;
}

