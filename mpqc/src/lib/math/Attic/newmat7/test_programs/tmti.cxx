
//#define WANT_STREAM

#include "include.h"

#include "newmatap.h"
//#include "newmatio.h"

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

void Clean(Matrix&, Real);

void WillNotConverge()
{
   Matrix A(10,10);
   Throw(ConvergenceException(A));
}

void trymati()
{
   Tracer et("Eighteenth test of Matrix package");
   Exception::PrintTrace(TRUE);
   ProgramException::SetAction(0);           // turn off error messages
   DataException::SetAction(0);
   ConvergenceException::SetAction(0);
   ColumnVector checks(14); checks = 1.0; checks(1) = 0.0;

   Try { WillNotConverge(); }
   Catch(ConvergenceException) { checks(2) = 0; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(1) = 1; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { Matrix M(10,10); SymmetricMatrix S = M; } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(3) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { Matrix M(10,10); M(10,11) = 2.0; } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(4) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { Matrix M(10,10); M = 0.0; M = M.i(); } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(1) = 1; }
   Catch(DataException) { checks(5) = 0; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { ColumnVector A(30), B(50);  A = 5; B = 3; FFT(A,B,A,B); } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(6) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try
   {
      ColumnVector A(30); A = 5; Matrix At = A.t();
      DiagonalMatrix D;
      SVD(At,D);
   } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(6) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { BandMatrix X(10,3,4); X(1,10) = 4.0; } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(7) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try
   {
      SymmetricMatrix S(10); S = 5;
      LowerTriangularMatrix L = Cholesky(S);
   } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(1) = 1; }
   Catch(DataException) { checks(8) = 0; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { BandMatrix M(10,3,5); M = 0.0; Matrix XM = M.i(); } 
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(1) = 1; }
   Catch(DataException) { checks(9) = 0; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;


   Try { ColumnVector X(10); ColumnVector Y; X = 5; X = X - Y; }
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(10) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try
   {
      UpperTriangularMatrix U(3); RowVector RV(3); RV = 10;
      U.Row(2) = RV;
   }
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(11) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { DiagonalMatrix D(3); D << 12 << 13 << 14 << 15; }
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(12) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;

   Try { ColumnVector D(3); D << 12 << 13; D << 1 << 2 << 3; }
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(13) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;


   Try {  { ColumnVector D(3); D << 12 << 13; }  }
   Catch(ConvergenceException) { checks(1) = 1; }
   Catch(InternalException) { checks(1) = 1; }
   Catch(ProgramException) { checks(14) = 0; }
   Catch(DataException) { checks(1) = 1; }
   Catch(SpaceException) { checks(1) = 1; }
   CatchAndThrow;


   ProgramException::SetAction(1);           // restore error messages
   DataException::SetAction(1);
   ConvergenceException::SetAction(1);

   Print(checks);

}
