
//#define WANT_STREAM

#include "include.h"
#include "newmatap.h"

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void Clean(Matrix&, Real);




void trymatd()
{
//   cout << "\nThirteenth test of Matrix package\n";
   Tracer et("Thirteenth test of Matrix package");
   Exception::PrintTrace(TRUE);
   Matrix X(5,20);
   int i,j;
   for (j=1;j<=20;j++) X(1,j) = j+1;
   for (i=2;i<=5;i++) for (j=1;j<=20; j++) X(i,j) = (long)X(i-1,j) * j % 1001;
   SymmetricMatrix S; S << X * X.t();
   Matrix SM = X * X.t() - S;
   Print(SM);
   LowerTriangularMatrix L = Cholesky(S);
   Matrix Diff = L*L.t()-S; Clean(Diff, 0.000000001);
   Print(Diff);
   LowerTriangularMatrix L1(5);
   HHDecompose(X,L1);
   Diff = L - L1; Clean(Diff,0.000000001); Print(Diff);

   ColumnVector C1(5);
   {
      Tracer et1("Stage 1");
      X.ReDimension(5,5);
      for (j=1;j<=5;j++) X(1,j) = j+1;
      for (i=2;i<=5;i++) for (j=1;j<=5; j++)
         X(i,j) = (long)X(i-1,j) * j % 1001;
      for (i=1;i<=5;i++) C1(i) = i*i;
      CroutMatrix A = X;
      ColumnVector C2 = A.i() * C1; C1 = X.i()  * C1;
      X = C1 - C2; Clean(X,0.000000001); Print(X); 
   }

   {
      Tracer et1("Stage 2");
      X.ReDimension(7,7);
      for (j=1;j<=7;j++) X(1,j) = j+1;
      for (i=2;i<=7;i++) for (j=1;j<=7; j++)
         X(i,j) = (long)X(i-1,j) * j % 1001;
      C1.ReDimension(7);
      for (i=1;i<=7;i++) C1(i) = i*i;
      RowVector R1 = C1.t();
      Diff = R1 * X.i() - ( X.t().i() * R1.t() ).t(); Clean(Diff,0.000000001);
      Print(Diff);
   }

   {
      Tracer et1("Stage 3");
      X.ReDimension(5,5);
      for (j=1;j<=5;j++) X(1,j) = j+1;
      for (i=2;i<=5;i++) for (j=1;j<=5; j++)
         X(i,j) = (long)X(i-1,j) * j % 1001;
      C1.ReDimension(5);
      for (i=1;i<=5;i++) C1(i) = i*i;
      CroutMatrix A1 = X*X;
      ColumnVector C2 = A1.i() * C1; C1 = X.i()  * C1; C1 = X.i()  * C1;
      X = C1 - C2; Clean(X,0.000000001); Print(X); 
   }


   {
      Tracer et1("Stage 4");
      int n = 40;
      SymmetricBandMatrix B(n,2); B = 0.0;
      for (i=1; i<=n; i++)
      {
	 B(i,i) = 6;
	 if (i<=n-1) B(i,i+1) = -4;
	 if (i<=n-2) B(i,i+2) = 1;
      }
      B(1,1) = 5; B(n,n) = 5;
      SymmetricMatrix A = B;
      ColumnVector X(n); X = 0; X(1) = 1;
      ColumnVector Y1 = A.i() * X;
      LowerTriangularMatrix C1 = Cholesky(A);
      ColumnVector Y2 = C1.t().i() * (C1.i() * X) - Y1;
      Clean(Y2, 0.000000001); Print(Y2);
  
      LowerBandMatrix C2 = Cholesky(B);
      Matrix M = C2 - C1; Clean(M, 0.000000001); Print(M);
      ColumnVector Y3 = C2.t().i() * (C2.i() * X) - Y1;
      Clean(Y3, 0.000000001); Print(Y3);
   }

   {
      Tracer et1("Stage 5");
      SymmetricMatrix A(4), B(4);
      A << 5
        << 1 << 4
        << 2 << 1 << 6
        << 1 << 0 << 1 << 7;
      B << 8
        << 1 << 5
        << 1 << 0 << 9
        << 2 << 1 << 0 << 6;
      LowerTriangularMatrix AB = Cholesky(A) * Cholesky(B);
      Matrix M = Cholesky(A) * B * Cholesky(A).t() - AB*AB.t();
      Clean(M, 0.000000001); Print(M);
      M = A * Cholesky(B); M = M * M.t() - A * B * A;
      Clean(M, 0.000000001); Print(M);
   }


//   cout << "\nEnd of Thirteenth test\n";
}
