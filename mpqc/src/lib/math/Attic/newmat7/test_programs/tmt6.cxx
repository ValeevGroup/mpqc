
//#define WANT_STREAM

#include "include.h"

#include "newmat.h"


/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void trymat6()
{
//   cout << "\nSixth test of Matrix package\n";
   Tracer et("Sixth test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i,j;


   DiagonalMatrix D(6);
   UpperTriangularMatrix U(6);
   for (i=1;i<=6;i++) { for (j=i;j<=6;j++) U(i,j)=i*i*i-50; D(i,i)=i*i+i-10; }
   LowerTriangularMatrix L=(U*3.0).t();
   SymmetricMatrix S(6);
   for (i=1;i<=6;i++) for (j=i;j<=6;j++) S(i,j)=i*i+2.0+j;
   Matrix MD=D; Matrix ML=L; Matrix MU=U; Matrix MS=S;
   Matrix M(6,6);
   for (i=1;i<=6;i++) for (j=1;j<=6;j++) M(i,j)=i*j+i*i-10.0;  
   {
      Tracer et1("Stage 1");
      Print(Matrix(MS+(-MS)));
      Print(Matrix((S+M)-(MS+M)));
      Print(Matrix((M+U)-(M+MU)));
      Print(Matrix((M+L)-(M+ML)));
   }
   {
      Tracer et1("Stage 2");
      Print(Matrix((M+D)-(M+MD)));
      Print(Matrix((U+D)-(MU+MD)));
      Print(Matrix((D+L)-(ML+MD)));
      Print(Matrix((-U+D)+MU-MD));
      Print(Matrix((-L+D)+ML-MD));
   }



//   cout << "\nEnd of sixth test\n";
}

