
//#define WANT_STREAM


#include "include.h"
#include "newmat.h"

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void Clean(Matrix&, Real);



void trymatc()
{
//   cout << "\nTwelfth test of Matrix package\n";
   Tracer et("Twelfth test of Matrix package");
   Exception::PrintTrace(TRUE);
   DiagonalMatrix D(15); D=1.5;
   Matrix A(15,15);
   int i,j;
   for (i=1;i<=15;i++) for (j=1;j<=15;j++) A(i,j)=i*i+j-150;
   { A = A + D; }
   ColumnVector B(15);
   for (i=1;i<=15;i++) B(i)=i+i*i-150.0;
   {
      Tracer et1("Stage 1");
      ColumnVector B1=B;
      B=(A*2.0).i() * B1;
      Matrix X = A*B-B1/2.0;
      Clean(X, 0.000000001); Print(X);
      A.ReDimension(3,5);
      for (i=1; i<=3; i++) for (j=1; j<=5; j++) A(i,j) = i+100*j;

      B = A.AsColumn()+10000;
      RowVector R = (A+10000).AsColumn().t();
      Print( RowVector(R-B.t()) );
   }

   {
      Tracer et1("Stage 2");
      B = A.AsColumn()+10000;
      Matrix XR = (A+10000).AsMatrix(15,1).t();
      Print( RowVector(XR-B.t()) );
   }

   {
      Tracer et1("Stage 3");
      B = (A.AsMatrix(15,1)+A.AsColumn())/2.0+10000;
      Matrix MR = (A+10000).AsColumn().t();
      Print( RowVector(MR-B.t()) );

      B = (A.AsMatrix(15,1)+A.AsColumn())/2.0;
      MR = A.AsColumn().t();
      Print( RowVector(MR-B.t()) );
   }

   {
      Tracer et1("Stage 4");
      B = (A.AsMatrix(15,1)+A.AsColumn())/2.0;
      RowVector R = A.AsColumn().t();
      Print( RowVector(R-B.t()) );
   }

   {
      Tracer et1("Stage 5");
      RowVector R = (A.AsColumn()-5000).t();
      B = ((R.t()+10000) - A.AsColumn())-5000;
      Print( RowVector(B.t()) );
   }

   {
      Tracer et1("Stage 6");
      B = A.AsColumn(); ColumnVector B1 = (A+10000).AsColumn() - 10000;
      Print(ColumnVector(B1-B));
   }

   {
      Tracer et1("Stage 7");
      Matrix X = B.AsMatrix(3,5); Print(Matrix(X-A));
      for (i=1; i<=3; i++) for (j=1; j<=5; j++) B(5*(i-1)+j) -= i+100*j;
      Print(B);
   }

   {
      Tracer et1("Stage 8");
      A.ReDimension(7,7); D.ReDimension(7);
      for (i=1; i<=7; i++) for (j=1; j<=7; j++) A(i,j) = i*j*j;
      for (i=1; i<=7; i++) D(i,i) = i;
      UpperTriangularMatrix U; U << A;
      Matrix X = A; for (i=1; i<=7; i++) X(i,i) = i;
      A.Inject(D); Print(Matrix(X-A));
      X = U; U.Inject(D); A = U; for (i=1; i<=7; i++) X(i,i) = i;
      Print(Matrix(X-A));
   }

   {
      Tracer et1("Stage 9");
      A.ReDimension(7,5);
      for (i=1; i<=7; i++) for (j=1; j<=5; j++) A(i,j) = i+100*j;
      Matrix Y = A; Y = Y - ((const Matrix&)A); Print(Y);
      Matrix X = A; // X.Release();
      Y = A; Y = ((const Matrix&)X) - A; Print(Y); Y = 0.0;
      Y = ((const Matrix&)X) - ((const Matrix&)A); Print(Y);
   }

//   cout << "\nEnd of twelfth test\n";
}
