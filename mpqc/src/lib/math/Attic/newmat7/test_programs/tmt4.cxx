
//#define WANT_STREAM


#include "include.h"

#include "newmat.h"


/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void trymat4()
{
//   cout << "\nFourth test of Matrix package\n";
   Tracer et("Fourth test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i,j;

   {
      Tracer et1("Stage 1");
      Matrix M(10,10);
      UpperTriangularMatrix U(10);
      for (i=1;i<=10;i++) for (j=1;j<=10;j++) M(i,j) = 100*i+j;
      U << -M;
      Matrix X1 = M.Rows(2,4);
      Matrix Y1 = U.t().Rows(2,4);
      Matrix X = U; { Print(Matrix(X.Columns(2,4).t()-Y1)); }
      RowVector RV = M.Row(5);
      {
         X.ReDimension(3,10);
         X.Row(1) << M.Row(2); X.Row(2) << M.Row(3); X.Row(3) << M.Row(4);
         Print(Matrix(X-X1));
      }
      {
         UpperTriangularMatrix V = U.SymSubMatrix(3,5);
         Matrix MV = U.SubMatrix(3,5,3,5); { Print(Matrix(MV-V)); }
         Matrix X2 = M.t().Columns(2,4); { Print(Matrix(X2-X1.t())); }
         Matrix Y2 = U.Columns(2,4); { Print(Matrix(Y2-Y1.t())); }
         ColumnVector CV = M.t().Column(5); { Print(ColumnVector(CV-RV.t())); }
         X.ReDimension(10,3); M = M.t();
         X.Column(1) << M.Column(2); X.Column(2) << M.Column(3);
         X.Column(3) << M.Column(4);
         Print(Matrix(X-X2));
      }
   }

   {
      Tracer et1("Stage 2");
      Matrix M; Matrix X; M.ReDimension(5,8);
      for (i=1;i<=5;i++) for (j=1;j<=8;j++) M(i,j) = 100*i+j;
      {
         X = M.Columns(5,8); M.Columns(5,8) << M.Columns(1,4);
             M.Columns(1,4) << X;
         X = M.Columns(3,4); M.Columns(3,4) << M.Columns(1,2);
             M.Columns(1,2) << X;
         X = M.Columns(7,8); M.Columns(7,8) << M.Columns(5,6);
             M.Columns(5,6) << X;
      }
      {
         X = M.Column(2); M.Column(2) = M.Column(1); M.Column(1) = X;
         X = M.Column(4); M.Column(4) = M.Column(3); M.Column(3) = X;
         X = M.Column(6); M.Column(6) = M.Column(5); M.Column(5) = X;
         X = M.Column(8); M.Column(8) = M.Column(7); M.Column(7) = X;
         X.ReDimension(5,8);
      }
      for (i=1;i<=5;i++) for (j=1;j<=8;j++) X(i,9-j) = 100*i+j;
      Print(Matrix(X-M));
   }

//   cout << "\nEnd of fourth test\n";
}

