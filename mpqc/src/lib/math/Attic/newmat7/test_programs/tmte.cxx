
//#define WANT_STREAM
#define WANT_MATH

#include "include.h"

#include "newmatap.h"

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void Clean(Matrix&, Real);
void Clean(DiagonalMatrix&, Real);




void trymate()
{
//   cout << "\nFourteenth test of Matrix package\n";
   Tracer et("Fourteenth test of Matrix package");
   Exception::PrintTrace(TRUE);

   {
      Tracer et1("Stage 1");
      Matrix A(8,5);
#ifndef ATandT
      Real   a[] =   { 22, 10,  2,  3,  7,
		       14,  7, 10,  0,  8,
		       -1, 13, -1,-11,  3,
		       -3, -2, 13, -2,  4,
			9,  8,  1, -2,  4,
			9,  1, -7,  5, -1,
			2, -6,  6,  5,  1,
			4,  5,  0, -2,  2 };
#else
      Real a[40];
      a[ 0]=22; a[ 1]=10; a[ 2]= 2; a[ 3]= 3; a[ 4]= 7;
      a[ 5]=14; a[ 6]= 7; a[ 7]=10; a[ 8]= 0; a[ 9]= 8;
      a[10]=-1; a[11]=13; a[12]=-1; a[13]=-11;a[14]= 3;
      a[15]=-3; a[16]=-2; a[17]=13; a[18]=-2; a[19]= 4;
      a[20]= 9; a[21]= 8; a[22]= 1; a[23]=-2; a[24]= 4;
      a[25]= 9; a[26]= 1; a[27]=-7; a[28]= 5; a[29]=-1;
      a[30]= 2; a[31]=-6; a[32]= 6; a[33]= 5; a[34]= 1;
      a[35]= 4; a[36]= 5; a[37]= 0; a[38]=-2; a[39]= 2;
#endif
      A << a;
      DiagonalMatrix D; Matrix U; Matrix V;
#ifdef ATandT
      int anc = A.Ncols(); DiagonalMatrix I(anc);     // AT&T 2.1 bug
#else
      DiagonalMatrix I(A.Ncols());
#endif
      I=1.0;
      SymmetricMatrix S1; S1 << A.t() * A;
      SymmetricMatrix S2; S2 << A * A.t();
      Real zero = 0.0; SVD(A+zero,D,U,V);
      DiagonalMatrix D1; SVD(A,D1); Print(DiagonalMatrix(D-D1));
      Matrix SU = U.t() * U - I; Clean(SU,0.000000001); Print(SU);
      Matrix SV = V.t() * V - I; Clean(SV,0.000000001); Print(SV);
      Matrix B = U * D * V.t() - A; Clean(B,0.000000001);Print(B);
      D1=0.0;  SVD(A,D1,A); Print(Matrix(A-U));
      SortDescending(D);
      D(1) -= sqrt(1248); D(2) -= 20; D(3) -= sqrt(384);
      Clean(D,0.000000001); Print(D);
      Jacobi(S1, D, V);
      V = S1 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      SortDescending(D); D(1)-=1248; D(2)-=400; D(3)-=384;
      Clean(D,0.000000001); Print(D);
      Jacobi(S2, D, V);
      V = S2 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      SortDescending(D); D(1)-=1248; D(2)-=400; D(3)-=384;
      Clean(D,0.000000001); Print(D);

      EigenValues(S1, D, V);
      V = S1 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      D(5)-=1248; D(4)-=400; D(3)-=384;
      Clean(D,0.000000001); Print(D);
      EigenValues(S2, D, V);
      V = S2 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      D(8)-=1248; D(7)-=400; D(6)-=384;
      Clean(D,0.000000001); Print(D);

      EigenValues(S1, D);
      D(5)-=1248; D(4)-=400; D(3)-=384;
      Clean(D,0.000000001); Print(D);
      EigenValues(S2, D);
      D(8)-=1248; D(7)-=400; D(6)-=384;
      Clean(D,0.000000001); Print(D);
   }
   {
      Tracer et1("Stage 2");
      Matrix A(20,21);
      for (int i=1; i<=20; i++) for (int j=1; j<=21; j++)
      { if (i>j) A(i,j) = 0; else if (i==j) A(i,j) = 21-i; else A(i,j) = -1; }
      A = A.t();
      SymmetricMatrix S1; S1 << A.t() * A;
      SymmetricMatrix S2; S2 << A * A.t();
      DiagonalMatrix D; Matrix U; Matrix V;
#ifdef ATandT
      int anc = A.Ncols(); DiagonalMatrix I(anc);     // AT&T 2.1 bug
#else
      DiagonalMatrix I(A.Ncols());
#endif
      I=1.0;
      SVD(A,D,U,V);
      Matrix SU = U.t() * U - I;    Clean(SU,0.000000001); Print(SU);
      Matrix SV = V.t() * V - I;    Clean(SV,0.000000001); Print(SV);
      Matrix B = U * D * V.t() - A; Clean(B,0.000000001);  Print(B);
      for (i=1; i<=20; i++)  D(i) -= sqrt((22-i)*(21-i));
      Clean(D,0.000000001); Print(D);
      Jacobi(S1, D, V);
      V = S1 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      SortDescending(D);
      for (i=1; i<=20; i++)  D(i) -= (22-i)*(21-i);
      Clean(D,0.000000001); Print(D);
      Jacobi(S2, D, V);
      V = S2 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      SortDescending(D);
      for (i=1; i<=20; i++)  D(i) -= (22-i)*(21-i);
      Clean(D,0.000000001); Print(D);

      EigenValues(S1, D, V);
      V = S1 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      for (i=1; i<=20; i++)  D(i) -= (i+1)*i;
      Clean(D,0.000000001); Print(D);
      EigenValues(S2, D, V);
      V = S2 - V * D * V.t(); Clean(V,0.000000001); Print(V);
      for (i=2; i<=21; i++)  D(i) -= (i-1)*i;
      Clean(D,0.000000001); Print(D);

      EigenValues(S1, D);
      for (i=1; i<=20; i++)  D(i) -= (i+1)*i;
      Clean(D,0.000000001); Print(D);
      EigenValues(S2, D);
      for (i=2; i<=21; i++)  D(i) -= (i-1)*i;
      Clean(D,0.000000001); Print(D);
   }

   {
      Tracer et1("Stage 3");
      Matrix A(30,30);
      for (int i=1; i<=30; i++) for (int j=1; j<=30; j++)
      { if (i>j) A(i,j) = 0; else if (i==j) A(i,j) = 1; else A(i,j) = -1; }
      Real d1 = A.LogDeterminant().Value();
      DiagonalMatrix D; Matrix U; Matrix V;
#ifdef ATandT
      int anc = A.Ncols(); DiagonalMatrix I(anc);     // AT&T 2.1 bug
#else
      DiagonalMatrix I(A.Ncols());
#endif
      I=1.0;
      SVD(A,D,U,V);
      Matrix SU = U.t() * U - I; Clean(SU,0.000000001); Print(SU);
      Matrix SV = V.t() * V - I; Clean(SV,0.000000001); Print(SV);
      Real d2 = D.LogDeterminant().Value();
      Matrix B = U * D * V.t() - A; Clean(B,0.000000001); Print(B);
      SortDescending(D);  // Print(D);
      Real d3 = D.LogDeterminant().Value();
      ColumnVector Test(3);
      Test(1) = d1 - 1; Test(2) = d2 - 1; Test(3) = d3 - 1;
      Clean(Test,0.00000001); Print(Test); // only 8 decimal figures
   }

//   cout << "\nEnd of Fourteenth test\n";
}
