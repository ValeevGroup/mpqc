
//#define WANT_STREAM


#include "include.h"

#include "newmat.h"


/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void trymat8()
{
//   cout << "\nEighth test of Matrix package\n";
   Tracer et("Eighth test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i;


   DiagonalMatrix D(6);
   for (i=1;i<=6;i++)  D(i,i)=i*i+i-10;
   DiagonalMatrix D2=D;
   Matrix MD=D;

   DiagonalMatrix D1(6); for (i=1;i<=6;i++) D1(i,i)=-100+i*i*i;
   Matrix MD1=D1;
   Print(Matrix(D*D1-MD*MD1));
   Print(Matrix((-D)*D1+MD*MD1));
   Print(Matrix(D*(-D1)+MD*MD1));
   DiagonalMatrix DX=D;
   {
      Tracer et1("Stage 1");
      DX=(DX+D1)*DX; Print(Matrix(DX-(MD+MD1)*MD));
      DX=D;
      DX=-DX*DX+(DX-(-D1))*((-D1)+DX);
#ifdef __GNUG__
      Matrix MX = Matrix(MD1);
      MD1=DX+(MX.t())*(MX.t()); Print(MD1);
#else
      MD1=DX+(Matrix(MD1).t())*(Matrix(MD1).t()); Print(MD1);
#endif
      DX=D; DX=DX; DX=D2-DX; Print(DiagonalMatrix(DX));
      DX=D;
   }
   {
      Tracer et1("Stage 2");
      D.Release(2);
      D1=D; D2=D;
      Print(DiagonalMatrix(D1-DX));
      Print(DiagonalMatrix(D2-DX));
      MD1=1.0;
      Print(Matrix(MD1-1.0));
   }

//   cout << "\nEnd of eighth test\n";
}
