
//#define WANT_STREAM


#include "include.h"

#include "newmat.h"


/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void trymat5()
{
//   cout << "\nFifth test of Matrix package\n";
   Tracer et("Fifth test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i,j;

   Matrix A(5,6);
   for (i=1;i<=5;i++) for (j=1;j<=6;j++) A(i,j)=1+i*j+i*i+j*j;
   ColumnVector CV(6);
   for (i=1;i<=6;i++) CV(i)=i*i+3;
   ColumnVector CV2(5); for (i=1;i<=5;i++) CV2(i)=1.0;
   ColumnVector CV1=CV;

   {
      CV=A*CV;
      RowVector RV=CV.t(); // RowVector RV; RV=CV.t();
      RV=RV-1.0;
      CV=(RV*A).t()+A.t()*CV2; CV1=(A.t()*A)*CV1 - CV;
      Print(CV1);
   }

   CV1.ReDimension(6);
   CV2.ReDimension(6);
   CV.ReDimension(6);
   for (i=1;i<=6;i++) { CV1(i)=i*3+1; CV2(i)=10-i; CV(i)=11+i*2; }
   ColumnVector CX=CV2-CV; { CX=CX+CV1; Print(CX); }
   Print(ColumnVector(CV1+CV2-CV));
   RowVector RV=CV.t(); RowVector RV1=CV1.t();
   RowVector R=RV-RV1; Print(RowVector(R-CV2.t()));

// test loading of list

   RV.ReDimension(10);
   for (i=1;i<=10;i++) RV(i) = i*i;
   RV1.ReDimension(10);
   RV1 << 1 << 4 << 9 << 16 << 25 << 36 << 49 << 64 << 81 << 100; // << 121;
   Print(RowVector(RV-RV1));
   et.ReName("Fifth test of Matrix package - at end");
   Matrix X(2,3);
   X << 11 << 12 << 13
     << 21 << 22 << 23;
   X(1,1) -= 11; X(1,2) -= 12; X(1,3) -= 13;
   X(2,1) -= 21; X(2,2) -= 22; X(2,3) -= 23;
   Print(X);
//   BandMatrix BM(7,2,2); BM << 23;

//   cout << "\nEnd of fifth test\n";
}
