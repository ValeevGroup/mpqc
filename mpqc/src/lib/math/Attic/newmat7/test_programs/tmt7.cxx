
//#define WANT_STREAM

#include "include.h"

#include "newmat.h"


/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void Clean(Matrix&, Real);

void trymat7()
{
//   cout << "\nSeventh test of Matrix package\n";
   Tracer et("Seventh test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i,j;


   DiagonalMatrix D(6);
   UpperTriangularMatrix U(6);
   for (i=1;i<=6;i++) { for (j=i;j<=6;j++) U(i,j)=i*i*j-50; D(i,i)=i*i+i-10; }
   LowerTriangularMatrix L=(U*3.0).t();
   SymmetricMatrix S(6);
   for (i=1;i<=6;i++) for (j=i;j<=6;j++) S(i,j)=i*i+2.0+j;
   Matrix MD=D; Matrix ML=L; Matrix MU=U;
   Matrix MS=S;
   Matrix M(6,6);
   for (i=1;i<=6;i++) for (j=1;j<=6;j++) M(i,j)=i*j+i*i-10.0;  
   {
      Tracer et1("Stage 1");
      Print(Matrix((S-M)-(MS-M)));
      Print(Matrix((-M-S)+(MS+M)));
      Print(Matrix((U-M)-(MU-M)));
   }
   {
      Tracer et1("Stage 2");
      Print(Matrix((L-M)+(M-ML)));
      Print(Matrix((D-M)+(M-MD)));
      Print(Matrix((D-S)+(MS-MD)));
      Print(Matrix((D-L)+(ML-MD)));
   }

   { M=MU.t(); }
   LowerTriangularMatrix LY=D.i()*U.t();
   {
      Tracer et1("Stage 3");
      MS=D*LY-M; Clean(MS,0.00000001); Print(MS);
      L=U.t();
      LY=D.i()*L; MS=D*LY-M; Clean(MS,0.00000001); Print(MS);
   }

//   cout << "\nEnd of seventh test\n";
}
