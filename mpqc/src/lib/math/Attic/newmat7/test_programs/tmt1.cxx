
//#define WANT_STREAM



#include "include.h"

#include "newmat.h"

/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);


void trymat1()
{
//   cout << "\nFirst test of Matrix package\n\n";
   Tracer et("First test of Matrix package");
   Exception::PrintTrace(TRUE);
   {
      Tracer et1("Stage 1");
      int i,j;

      LowerTriangularMatrix L(10);
      for (i=1;i<=10;i++) for (j=1;j<=i;j++) L(i,j)=2.0+i*i+j;
      SymmetricMatrix S(10);
      for (i=1;i<=10;i++) for (j=1;j<=i;j++) S(i,j)=i*j+1.0;
      SymmetricMatrix S1 = S / 2.0;
      S = S1 * 2.0;
      UpperTriangularMatrix U=L.t()*2.0;
      Print(LowerTriangularMatrix(L-U.t()*0.5));
      DiagonalMatrix D(10);
      for (i=1;i<=10;i++) D(i,i)=(i-4)*(i-5)*(i-6);
      Matrix M=(S+U-D+L)*(L+U-D+S);
      DiagonalMatrix DD=D*D;
      LowerTriangularMatrix LD=L*D;
      // expressions split for Turbo C
      Matrix M1 = S*L + U*L - D*L + L*L + 10.0;
      { M1 = M1 + S*U + U*U - D*U + L*U - S*D; }
      { M1 = M1 - U*D + DD - LD + S*S; }
      { M1 = M1 + U*S - D*S + L*S - 10.0; }
      M=M1-M;
      Print(M);
   }
   {
      Tracer et1("Stage 2");
      int i,j;

      LowerTriangularMatrix L(9);
      for (i=1;i<=9;i++) for (j=1;j<=i;j++) L(i,j)=1.0+j;
      UpperTriangularMatrix U1(9);
      for (j=1;j<=9;j++) for (i=1;i<=j;i++) U1(i,j)=1.0+i;
      LowerTriangularMatrix LX(9);
      for (i=1;i<=9;i++) for (j=1;j<=i;j++) LX(i,j)=1.0+i*i;
      UpperTriangularMatrix UX(9);
      for (j=1;j<=9;j++) for (i=1;i<=j;i++) UX(i,j)=1.0+j*j;
      {
         L=L+LX/0.5;   L=L-LX*3.0;   L=LX*2.0+L;
         U1=U1+UX*2.0; U1=U1-UX*3.0; U1=UX*2.0+U1;
      }


      SymmetricMatrix S(9);
      for (i=1;i<=9;i++) for (j=1;j<=i;j++) S(i,j)=i*i+j;
      {
         SymmetricMatrix S1 = S;
         S=S1+5.0;
         S=S-3.0;
      }

      DiagonalMatrix D(9);
      for (i=1;i<=9;i++) D(i,i)=S(i,i);
      UpperTriangularMatrix U=L.t()*2.0;
      {
         U1=U1*2.0 - U;  Print(U1);
         L=L*2.0-D; U=U-D;
      }
      Matrix M=U+L; S=S*2.0; M=S-M; Print(M);
   }
//   cout << "\nEnd of first test\n";
}

