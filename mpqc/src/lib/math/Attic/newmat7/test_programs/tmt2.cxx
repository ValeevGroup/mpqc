
//#define WANT_STREAM


#include "include.h"

#include "newmat.h"

/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);


void trymat2()
{
//   cout << "\nSecond test of Matrix package\n\n";
   Tracer et("Second test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i,j;
   Matrix M(3,5);
   for (i=1; i<=3; i++) for (j=1; j<=5; j++) M(i,j) = 100*i + j;
   Matrix X(8,10);
   for (i=1; i<=8; i++) for (j=1; j<=10; j++) X(i,j) = 1000*i + 10*j;
   Matrix Y = X; Matrix Z = X;
   { X.SubMatrix(2,4,3,7) << M; }
   for (i=1; i<=3; i++) for (j=1; j<=5; j++) Y(i+1,j+2) = 100*i + j;
   Print(Matrix(X-Y));


   Real a[15]; Real* r = a;
   for (i=1; i<=3; i++) for (j=1; j<=5; j++) *r++ = 100*i + j;
   { Z.SubMatrix(2,4,3,7) << a; }
   Print(Matrix(Z-Y));

   { M=33; X.SubMatrix(2,4,3,7) << M; }
   { Z.SubMatrix(2,4,3,7) = 33; }
   Print(Matrix(Z-X));

   for (i=1; i<=8; i++) for (j=1; j<=10; j++) X(i,j) = 1000*i + 10*j;
   Y = X;
   UpperTriangularMatrix U(5);
   for (i=1; i<=5; i++) for (j=i; j<=5; j++) U(i,j) = 100*i + j;
   { X.SubMatrix(3,7,5,9) << U; }
   for (i=1; i<=5; i++) for (j=i; j<=5; j++) Y(i+2,j+4) = 100*i + j;
   for (i=1; i<=5; i++) for (j=1; j<i; j++) Y(i+2,j+4) = 0.0;
   Print(Matrix(X-Y));
   for (i=1; i<=8; i++) for (j=1; j<=10; j++) X(i,j) = 1000*i + 10*j;
   Y = X;
   for (i=1; i<=5; i++) for (j=i; j<=5; j++) U(i,j) = 100*i + j;
   { X.SubMatrix(3,7,5,9).Inject(U); }
   for (i=1; i<=5; i++) for (j=i; j<=5; j++) Y(i+2,j+4) = 100*i + j;
   Print(Matrix(X-Y));


   // test growing and shrinking a vector
   ColumnVector V(100);
   for (i=1;i<=100;i++) V(i) = i*i+i;
   V = V.Rows(1,50);               // to get first 50 vlaues.

   {
      V.Release(); ColumnVector VX=V;
      V.ReDimension(100); V = 0.0; V.Rows(1,50)=VX;
   }                               // V now length 100

#if __ZTC__
   ColumnVector VM = V; VM = 100;
#else
   M=V; M=100;                     // to make sure V will hold its values
#endif
   for (i=1;i<=50;i++) V(i) -= i*i+i;
   Print(V);


//   cout << "\nEnd of second test\n";
}
