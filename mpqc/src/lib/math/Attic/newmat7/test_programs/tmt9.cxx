
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
void Clean(DiagonalMatrix&, Real);

void trymat9()
{
//   cout << "\nNinth test of Matrix package\n";
   Tracer et("Ninth test of Matrix package");
   Exception::PrintTrace(TRUE);


   int i; int j;
   Matrix A(7,7); Matrix X(7,3);
   for (i=1;i<=7;i++) for (j=1;j<=7;j++) A(i,j)=i*i+j+((i==j) ? 1 : 0);
   for (i=1;i<=7;i++) for (j=1;j<=3;j++) X(i,j)=i-j;
   Matrix B = A.i(); DiagonalMatrix D(7); D=1.0;
   {
      Tracer et1("Stage 1");
      Matrix Q = B*A-D; Clean(Q, 0.000000001); Print(Q);
      Q=A; Q = Q.i() * X; Q = A*Q - X; Clean(Q, 0.000000001); Print(Q);
      Q=X; Q = A.i() * Q; Q = A*Q - X; Clean(Q, 0.000000001); Print(Q);
   }
   for (i=1;i<=7;i++) D(i,i)=i*i+1;
   DiagonalMatrix E(3); for (i=1;i<=3;i++) E(i,i)=i+23;
   {
      Tracer et1("Stage 2");
      Matrix DXE = D.i() * X * E;
      DXE = E.i() * DXE.t() * D - X.t(); Clean(DXE, 0.00000001); Print(DXE); 
      E=D; for (i=1;i<=7;i++) E(i,i)=i*3+1;
   }
   DiagonalMatrix F=D;
   {
      Tracer et1("Stage 3");
      F=E.i()*F; F=F*E-D; Clean(F,0.00000001); Print(F);
      F=E.i()*D; F=F*E-D; Clean(F,0.00000001); Print(F);
   }
   { F=E; F=F.i()*D; F=F*E-D; Clean(F,0.00000001); Print(F); }
//   cout << "\nEnd of ninth test\n";
}

