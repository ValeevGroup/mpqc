//$$ example.cxx                             Example of use of matrix package

#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns

#include "include.h"                 // include standard files

#include "newmatap.h"                // need matrix applications

Real t3(Real);                       // round to 3 decimal places


// demonstration of matrix package on linear regression problem


void test1(Real* y, Real* x1, Real* x2, int nobs, int npred)
{
   cout << "\n\nTest 1 - traditional\n";

   // traditional sum of squares and products method of calculation
   // with subtraction of means

   // make matrix of predictor values
   Matrix X(nobs,npred);

   // load x1 and x2 into X
   //    [use << rather than = with submatrices and/or loading arrays]
   X.Column(1) << x1;  X.Column(2) << x2;

   // vector of Y values
   ColumnVector Y(nobs); Y << y;

   // make vector of 1s
   ColumnVector Ones(nobs); Ones = 1.0;

   // calculate means (averages) of x1 and x2 [ .t() takes transpose]
   RowVector M = Ones.t() * X / nobs;

   // and subtract means from x1 and x1
   Matrix XC(nobs,npred);
   XC = X - Ones * M;

   // do the same to Y [need "Real" to convert 1x1 matrix to scalar]
   ColumnVector YC(nobs);
   Real m = (Ones.t() * Y).AsScalar() / nobs;  YC = Y - Ones * m;

   // form sum of squares and product matrix
   //    [use << rather than = for copying Matrix into SymmetricMatrix]
   SymmetricMatrix SSQ; SSQ << XC.t() * XC;

   // calculate estimate
   //    [bracket last two terms to force this multiplication first]
   //    [ .i() means inverse, but inverse is not explicity calculated]
   ColumnVector A = SSQ.i() * (XC.t() * YC);

   // calculate estimate of constant term
   Real a = m - (M * A).AsScalar();

   // Get variances of estimates from diagonal elements of invoice of SSQ
   //    [ we are taking inverse of SSQ; would have been better to use
   //        CroutMatrix method - see documentation ]
   Matrix ISSQ = SSQ.i(); DiagonalMatrix D; D << ISSQ;
   ColumnVector V = D.AsColumn();
   Real v = 1.0/nobs + (M * ISSQ * M.t()).AsScalar();
					    // for calc variance const

   // Calculate fitted values and residuals
   int npred1 = npred+1;
   ColumnVector Fitted = X * A + a;
   ColumnVector Residual = Y - Fitted;
   Real ResVar = Residual.SumSquare() / (nobs-npred1);

   // Get diagonals of Hat matrix (an expensive way of doing this)
   Matrix X1(nobs,npred1); X1.Column(1)<<Ones; X1.Columns(2,npred1)<<X;
   DiagonalMatrix Hat;  Hat << X1 * (X1.t() * X1).i() * X1.t();

   // print out answers
   cout << "\nEstimates and their standard errors\n\n";
   cout << a <<"\t"<< sqrt(v*ResVar) << "\n";
   for (int i=1; i<=npred; i++)
   cout << A(i) <<"\t"<< sqrt(V(i)*ResVar) << "\n";
   cout << "\nObservations, fitted value, residual value, hat value\n";
   for (i=1; i<=nobs; i++)
      cout << X(i,1) <<"\t"<< X(i,2) <<"\t"<< Y(i) <<"\t"<<
      t3(Fitted(i)) <<"\t"<< t3(Residual(i)) <<"\t"<< t3(Hat(i)) <<"\n";
   cout << "\n\n";
}

void test2(Real* y, Real* x1, Real* x2, int nobs, int npred)
{
   cout << "\n\nTest 2 - Cholesky\n";

   // traditional sum of squares and products method of calculation
   // with subtraction of means - using Cholesky decomposition

   Matrix X(nobs,npred);
   X.Column(1) << x1;  X.Column(2) << x2;
   ColumnVector Y(nobs); Y << y;
   ColumnVector Ones(nobs); Ones = 1.0;
   RowVector M = Ones.t() * X / nobs;
   Matrix XC(nobs,npred);
   XC = X - Ones * M;
   ColumnVector YC(nobs);
   Real m = (Ones.t() * Y).AsScalar() / nobs;  YC = Y - Ones * m;
   SymmetricMatrix SSQ; SSQ << XC.t() * XC;

   // Cholesky decomposition of SSQ
   LowerTriangularMatrix L = Cholesky(SSQ);

   // calculate estimate
   ColumnVector A = L.t().i() * (L.i() * (XC.t() * YC));

   // calculate estimate of constant term
   Real a = m - (M * A).AsScalar();

   // Get variances of estimates from diagonal elements of invoice of SSQ
   DiagonalMatrix D; D << L.t().i() * L.i();
   ColumnVector V = D.AsColumn();
   Real v = 1.0/nobs + (L.i() * M.t()).SumSquare();

   // Calculate fitted values and residuals
   int npred1 = npred+1;
   ColumnVector Fitted = X * A + a;
   ColumnVector Residual = Y - Fitted;
   Real ResVar = Residual.SumSquare() / (nobs-npred1);

   // Get diagonals of Hat matrix (an expensive way of doing this)
   Matrix X1(nobs,npred1); X1.Column(1)<<Ones; X1.Columns(2,npred1)<<X;
   DiagonalMatrix Hat;  Hat << X1 * (X1.t() * X1).i() * X1.t();

   // print out answers
   cout << "\nEstimates and their standard errors\n\n";
   cout << a <<"\t"<< sqrt(v*ResVar) << "\n";
   for (int i=1; i<=npred; i++)
      cout << A(i) <<"\t"<< sqrt(V(i)*ResVar) << "\n";
   cout << "\nObservations, fitted value, residual value, hat value\n";
   for (i=1; i<=nobs; i++)
      cout << X(i,1) <<"\t"<< X(i,2) <<"\t"<< Y(i) <<"\t"<<
      t3(Fitted(i)) <<"\t"<< t3(Residual(i)) <<"\t"<< t3(Hat(i)) <<"\n";
   cout << "\n\n";
}

void test3(Real* y, Real* x1, Real* x2, int nobs, int npred)
{
   cout << "\n\nTest 3 - Householder triangularisation\n";

   // Householder triangularisation method
 
   // load data - 1s into col 1 of matrix
   int npred1 = npred+1;
   Matrix X(nobs,npred1); ColumnVector Y(nobs);
   X.Column(1) = 1.0;  X.Column(2) << x1;  X.Column(3) << x2;  Y << y;

   // do Householder triangularisation
   // no need to deal with constant term separately
   Matrix XT = X.t();             // Want data to be along rows
   RowVector YT = Y.t();
   LowerTriangularMatrix L; RowVector M;
   HHDecompose(XT, L); HHDecompose(XT, YT, M); // YT now contains resids
   ColumnVector A = L.t().i() * M.t();
   ColumnVector Fitted = X * A;
   Real ResVar = YT.SumSquare() / (nobs-npred1);

   // get variances of estimates
   L = L.i(); DiagonalMatrix D; D << L.t() * L;

   // Get diagonals of Hat matrix
   DiagonalMatrix Hat;  Hat << XT.t() * XT;

   // print out answers
   cout << "\nEstimates and their standard errors\n\n";
   for (int i=1; i<=npred1; i++)
      cout << A(i) <<"\t"<< sqrt(D(i)*ResVar) << "\n";
   cout << "\nObservations, fitted value, residual value, hat value\n";
   for (i=1; i<=nobs; i++)
      cout << X(i,2) <<"\t"<< X(i,3) <<"\t"<< Y(i) <<"\t"<<
      t3(Fitted(i)) <<"\t"<< t3(YT(i)) <<"\t"<< t3(Hat(i)) <<"\n";
   cout << "\n\n";
}

void test4(Real* y, Real* x1, Real* x2, int nobs, int npred)
{
   cout << "\n\nTest 4 - singular value\n";

   // Singular value decomposition method
 
   // load data - 1s into col 1 of matrix
   int npred1 = npred+1;
   Matrix X(nobs,npred1); ColumnVector Y(nobs);
   X.Column(1) = 1.0;  X.Column(2) << x1;  X.Column(3) << x2;  Y << y;

   // do SVD
   Matrix U, V; DiagonalMatrix D;
   SVD(X,D,U,V);                              // X = U * D * V.t()
   ColumnVector Fitted = U.t() * Y;
   ColumnVector A = V * ( D.i() * Fitted );
   Fitted = U * Fitted;
   ColumnVector Residual = Y - Fitted;
   Real ResVar = Residual.SumSquare() / (nobs-npred1);

   // get variances of estimates
   D << V * (D * D).i() * V.t();

   // Get diagonals of Hat matrix
   DiagonalMatrix Hat;  Hat << U * U.t();

   // print out answers
   cout << "\nEstimates and their standard errors\n\n";
   for (int i=1; i<=npred1; i++)
      cout << A(i) <<"\t"<< sqrt(D(i)*ResVar) << "\n";
   cout << "\nObservations, fitted value, residual value, hat value\n";
   for (i=1; i<=nobs; i++)
      cout << X(i,2) <<"\t"<< X(i,3) <<"\t"<< Y(i) <<"\t"<<
      t3(Fitted(i)) <<"\t"<< t3(Residual(i)) <<"\t"<< t3(Hat(i)) <<"\n";
   cout << "\n\n";
}

main()
{
   cout << "\nDemonstration of Matrix package\n\n";

   // Test for any memory not deallocated after running this program
   Real* s1; { ColumnVector A(8000); s1 = A.Store(); }

   {
      // the data

#ifndef ATandT
      Real y[]  = { 8.3, 5.5, 8.0, 8.5, 5.7, 4.4, 6.3, 7.9, 9.1 };
      Real x1[] = { 2.4, 1.8, 2.4, 3.0, 2.0, 1.2, 2.0, 2.7, 3.6 };
      Real x2[] = { 1.7, 0.9, 1.6, 1.9, 0.5, 0.6, 1.1, 1.0, 0.5 };
#else             // for compilers that don't understand aggregrates
      Real y[9], x1[9], x2[9];
      y[0]=8.3; y[1]=5.5; y[2]=8.0; y[3]=8.5; y[4]=5.7;
      y[5]=4.4; y[6]=6.3; y[7]=7.9; y[8]=9.1;
      x1[0]=2.4; x1[1]=1.8; x1[2]=2.4; x1[3]=3.0; x1[4]=2.0;
      x1[5]=1.2; x1[6]=2.0; x1[7]=2.7; x1[8]=3.6;
      x2[0]=1.7; x2[1]=0.9; x2[2]=1.6; x2[3]=1.9; x2[4]=0.5;
      x2[5]=0.6; x2[6]=1.1; x2[7]=1.0; x2[8]=0.5;
#endif

      int nobs = 9;                           // number of observations
      int npred = 2;                          // number of predictor values

      // we want to find the values of a,a1,a2 to give the best
      // fit of y[i] with a0 + a1*x1[i] + a2*x2[i]
      // Also print diagonal elements of hat matrix, X*(X.t()*X).i()*X.t()

      // this example demonstrates four methods of calculation

      Try
      {
         test1(y, x1, x2, nobs, npred);
         test2(y, x1, x2, nobs, npred);
         test3(y, x1, x2, nobs, npred);
         test4(y, x1, x2, nobs, npred);
      }
      Catch(DataException) { cout << "\nInvalid data\n"; }
      Catch(SpaceException) { cout << "\nMemory exhausted\n"; }
      CatchAll { cout << "\nUnexpected program failure\n"; }
   }

#ifdef DO_FREE_CHECK
   FreeCheck::Status();
#endif
   Real* s2; { ColumnVector A(8000); s2 = A.Store(); }
   cout << "\n\nChecking for lost memory: "
      << (unsigned long)s1 << " " << (unsigned long)s2 << " ";
   if (s1 != s2) cout << " - error\n"; else cout << " - ok\n";

   return 0;

}

Real t3(Real r) { return int(r*1000) / 1000.0; }
