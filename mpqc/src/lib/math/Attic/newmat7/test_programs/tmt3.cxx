
//#define WANT_STREAM

#include "include.h"

#include "newmat.h"


/**************************** test program ******************************/

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void trymat3()
{
//   cout << "\nThird test of Matrix package\n";
   Tracer et("Third test of Matrix package");
   Exception::PrintTrace(TRUE);

   int i,j;

   {
      Tracer et1("Stage 1");
      SymmetricMatrix S(7);
      for (i=1;i<=7;i++) for (j=1;j<=i;j++) S(i,j)=i*i+j;
      S=-S+2.0;

      DiagonalMatrix D(7);
      for (i=1;i<=7;i++) D(i,i)=S(i,i);

      Matrix M4(7,7); { M4=D+(D+4.0); M4=M4-D*2.0;  M4=M4-4.0; Print(M4); }
      SymmetricMatrix S2=D; Matrix M2=S2;  { M2=-D+M2; Print(M2); }
      UpperTriangularMatrix U2=D; { M2=U2; M2=D-M2; Print(M2); }
      LowerTriangularMatrix L2=D; { M2=L2; M2=D-M2; Print(M2); }
      M2=D; M2=M2-D; Print(M2);
      for (i=1;i<=7;i++) for (j=1;j<=i;j++) L2(i,j)=2.0-i*i-j;
      U2=L2.t(); D=D.t(); S=S.t();
      M4=(L2-1.0)+(U2+1.0)-D-S; Print(M4);
      M4=(-L2+1.0)+(-U2-1.0)+D+S; Print(M4);
   }

   {
      Tracer et1("Stage 2");
      DiagonalMatrix D(6);
      for (i=1;i<=6;i++) D(i,i)=i*3.0+i*i+2.0;
      UpperTriangularMatrix U2(7); LowerTriangularMatrix L2(7);
      for (i=1;i<=7;i++) for (j=1;j<=i;j++) L2(i,j)=2.0-i*i+j;
      { U2=L2.t(); }
      DiagonalMatrix D1(7); for (i=1;i<=7;i++) D1(i,i)=(i-2)*(i-4);
      Matrix M2(6,7);
      for (i=1;i<=6;i++) for (j=1;j<=7;j++) M2(i,j)=2.0+i*j+i*i+2.0*j*j;
      Matrix MD=D; SymmetricMatrix MD1(1); MD1=D1;
      Matrix MX=MD*M2*MD1 - D*(M2*D1); Print(MX);
      MX=MD*M2*MD1 - (D*M2)*D1; Print(MX);
      {
         D.ReDimension(7); for (i=1;i<=7;i++) D(i,i)=i*3.0+i*i+2.0;
         LowerTriangularMatrix LD(1); LD=D;
         UpperTriangularMatrix UD(1); UD=D;
         M2=U2; M2=LD*M2*MD1 - D*(U2*D1); Print(M2);
         M2=U2; M2=UD*M2*MD1 - (D*U2)*D1; Print(M2);
         M2=L2; M2=LD*M2*MD1 - D*(L2*D1); Print(M2);
         M2=L2; M2=UD*M2*MD1 - (D*L2)*D1; Print(M2);
      }
   }

//   cout << "\nEnd of third test\n";
}

