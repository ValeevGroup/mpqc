
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


static void SlowFT(const ColumnVector& a, const ColumnVector&b,
   ColumnVector& x, ColumnVector& y)
{
   int n = a.Nrows();
   x.ReDimension(n); y.ReDimension(n);
   Real f = 6.2831853071795864769/n;
   for (int j=1; j<=n; j++)
   {
      Real sumx = 0.0; Real sumy = 0.0;
      for (int k=1; k<=n; k++)
      {
	 Real theta = - (j-1) * (k-1) * f;
	 Real c = cos(theta);  Real s = sin(theta);
	 sumx += c * a(k) - s * b(k);
	 sumy += s * a(k) + c * b(k);
      }
      x(j) = sumx; y(j) = sumy;
   }
}

static void test(int n)
{
   Tracer et("Test FFT");

   ColumnVector A(n), B(n), X, Y;
   for (int i=0; i<n; i++)
   {
      Real f = 6.2831853071795864769*i/n;
      A.element(i) = fabs(sin(7.0*f) + 0.5 * cos(19.0 * f)) + (Real)i/(Real)n;
      B.element(i) = fabs(0.25 * cos(31.0 * f)) + (Real)i/(Real)n;
   }
   FFT(A, B, X, Y);
   FFTI(X, Y, X, Y);
   X = X - A; Y = Y - B;
   Clean(X,0.000000001); Clean(Y,0.000000001); Print(X); Print(Y);
}

static void test1(int n)
{
   Tracer et("Test RealFFT");

   ColumnVector A(n), B(n), X, Y;
   for (int i=0; i<n; i++)
   {
      Real f = 6.2831853071795864769*i/n;
      A.element(i) = fabs(sin(7.0*f) + 0.5 * cos(19.0 * f)) + (Real)i/(Real)n;
   }
   B = 0.0;
   FFT(A, B, X, Y);
   int n2 = n/2+1;
   ColumnVector X1,Y1,X2,Y2;
   RealFFT(A, X1, Y1);
   X2 = X1 - X.Rows(1,n2); Y2 = Y1 - Y.Rows(1,n2);
   Clean(X2,0.000000001); Clean(Y2,0.000000001); Print(X2); Print(Y2);
   RealFFTI(X1,Y1,B);
   B = A - B;
   Clean(B,0.000000001); Print(B);
}

static void test2(int n)
{
   Tracer et("cf FFT and slow FT");

   ColumnVector A(n), B(n), X, Y, X1, Y1;
   for (int i=0; i<n; i++)
   {
      Real f = 6.2831853071795864769*i/n;
      A.element(i) = fabs(sin(7.0*f) - 0.5 * cos(19.0 * f)) + (Real)i/(Real)n;
      B.element(i) = fabs(0.25 * cos(31.0 * f)) - (Real)i/(Real)n;
   }
   FFT(A, B, X, Y);
   SlowFT(A, B, X1, Y1);
   X = X - X1; Y = Y - Y1;
   Clean(X,0.000000001); Clean(Y,0.000000001); Print(X); Print(Y);
}

void trymatf()
{
   Tracer et("Fifteenth test of Matrix package");
   Exception::PrintTrace(TRUE);
//   cout << "\nFifteenth test of Matrix package\n";
//   cout << "\n";

   ColumnVector A(12), B(12);
   for (int i = 1; i <=12; i++)
   {
      Real i1 = i - 1;
      A(i) = .7
	   + .2 * cos(6.2831853071795864769 * 4.0 * i1 / 12)
	   + .3 * sin(6.2831853071795864769 * 3.0 * i1 / 12);
      B(i) = .9
	   + .5 * sin(6.2831853071795864769 * 2.0 * i1 / 12)
	   + .4 * cos(6.2831853071795864769 * 1.0 * i1 / 12);
   }
   FFT(A, B, A, B);
   A(1) -= 8.4; A(3) -= 3.0; A(5) -= 1.2; A(9) -= 1.2; A(11) += 3.0;
   B(1) -= 10.8; B(2) -= 2.4; B(4) += 1.8; B(10) -= 1.8; B(12) -= 2.4;
   Clean(A,0.000000001); Clean(B,0.000000001); Print(A); Print(B);



   test(2048);
   test(2000);
   test(27*81);
   test(2310);
   test(49*49);
   test(13*13*13);
   test(43*47);
   test1(98);
   test1(100);
   test1(2048);
   test1(2000);
   test1(35*35*2);
   test2(13);
   test2(12);
   test2(9);
   test2(16);

//   cout << "\nEnd of Fifteenth test\n";
}
