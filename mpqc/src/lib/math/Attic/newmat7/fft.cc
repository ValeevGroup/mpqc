//$$ fft.cxx                         Fast fourier transform

// Copyright (C) 1991,2,3: R B Davies


#define WANT_MATH

#include "include.h"

#include "newmatap.h"


static void cossin(int n, int d, Real& c, Real& s)
// calculate cos(twopi*n/d) and sin(twopi*n/d)
// minimise roundoff error
{
   long n4 = n * 4; int sector = (int)floor( (Real)n4 / (Real)d + 0.5 );
   n4 -= sector * d;
   if (sector < 0) sector = 3 - (3 - sector) % 4; else sector %= 4;
   Real ratio = 1.5707963267948966192 * (Real)n4 / (Real)d;

   switch (sector)
   {
   case 0: c =  cos(ratio); s =  sin(ratio); break;
   case 1: c = -sin(ratio); s =  cos(ratio); break;
   case 2: c = -cos(ratio); s = -sin(ratio); break;
   case 3: c =  sin(ratio); s = -cos(ratio); break;
   }
}

static void fftstep(ColumnVector& A, ColumnVector& B, ColumnVector& X,
   ColumnVector& Y, int after, int now, int before)
{
   Tracer trace("FFT(step)");
   // const Real twopi = 6.2831853071795864769;
   const int gamma = after * before;  const int delta = now * after;
   // const Real angle = twopi / delta;  Real temp;
   // Real r_omega = cos(angle);  Real i_omega = -sin(angle);
   Real r_arg = 1.0;  Real i_arg = 0.0;
   Real* x = X.Store();  Real* y = Y.Store();   // pointers to array storage
   const int m = A.Nrows() - gamma;

   for (int j = 0; j < now; j++)
   {
      Real* a = A.Store(); Real* b = B.Store(); // pointers to array storage
      Real* x1 = x; Real* y1 = y; x += after; y += after;
      for (int ia = 0; ia < after; ia++)
      {
	 // generate sins & cosines explicitly rather than iteratively
	 // for more accuracy; but slower
	 cossin(-(j*after+ia), delta, r_arg, i_arg);

	 Real* a1 = a++; Real* b1 = b++; Real* x2 = x1++; Real* y2 = y1++;
	 if (now==2)
	 {
	    int ib = before; while (ib--)
	    {
	       Real* a2 = m + a1; Real* b2 = m + b1; a1 += after; b1 += after;
	       Real r_value = *a2; Real i_value = *b2;
	       *x2 = r_value * r_arg - i_value * i_arg + *(a2-gamma);
	       *y2 = r_value * i_arg + i_value * r_arg + *(b2-gamma);
	       x2 += delta; y2 += delta;
	    }
	 }
	 else
	 {
	    int ib = before; while (ib--)
	    {
	       Real* a2 = m + a1; Real* b2 = m + b1; a1 += after; b1 += after;
	       Real r_value = *a2; Real i_value = *b2;
	       int in = now-1; while (in--)
	       {
		  // it should be possible to make this faster
		  // hand code for now = 2,3,4,5,8
		  // use symmetry to halve number of operations
		  a2 -= gamma; b2 -= gamma;  Real temp = r_value;
		  r_value = r_value * r_arg - i_value * i_arg + *a2;
		  i_value = temp    * i_arg + i_value * r_arg + *b2;
	       }
	       *x2 = r_value; *y2 = i_value;   x2 += delta; y2 += delta;
	    }
	 }

         // temp = r_arg;
         // r_arg = r_arg * r_omega - i_arg * i_omega;
         // i_arg = temp  * i_omega + i_arg * r_omega;

      }
   }
}


void FFT(const ColumnVector& U, const ColumnVector& V,
   ColumnVector& X, ColumnVector& Y)
{
   // from Carl de Boor (1980), Siam J Sci Stat Comput, 1 173-8
   Tracer trace("FFT");
   const int n = U.Nrows();                     // length of arrays
   if (n != V.Nrows() || n == 0)
      Throw(ProgramException("Vector lengths unequal or zero", U, V));
   ColumnVector A = U; ColumnVector B = V;
   X.ReDimension(n); Y.ReDimension(n);
   const int nextmx = 8;
#ifndef ATandT
   int prime[8] = { 2,3,5,7,11,13,17,19 };
#else
   int prime[8];
   prime[0]=2; prime[1]=3; prime[2]=5; prime[3]=7;
   prime[4]=11; prime[5]=13; prime[6]=17; prime[7]=19;
#endif
   int after = 1; int before = n; int next = 0; Boolean inzee = TRUE;

   do
   {
      int now, b1;
      for (;;)
      {
	 if (next < nextmx) now = prime[next];
	 b1 = before / now;  if (b1 * now == before) break;
	 next++; now += 2;
      }
      before = b1;

      if (inzee) fftstep(A, B, X, Y, after, now, before);
      else fftstep(X, Y, A, B, after, now, before);

      inzee = !inzee; after *= now;
   }
   while (before != 1);

   if (inzee) { A.Release(); X = A; B.Release(); Y = B; }
}


void FFTI(const ColumnVector& U, const ColumnVector& V,
   ColumnVector& X, ColumnVector& Y)
{
   // Inverse transform
   Tracer trace("FFTI");
   FFT(U,-V,X,Y);
   const Real n = X.Nrows(); X = X / n; Y = Y / (-n);
}

void RealFFT(const ColumnVector& U, ColumnVector& X, ColumnVector& Y)
{
   // Fourier transform of a real series
   Tracer trace("RealFFT");
   const int n = U.Nrows();                     // length of arrays
   const int n2 = n / 2;
   if (n != 2 * n2)
      Throw(ProgramException("Vector length not multiple of 2", U));
   ColumnVector A(n2), B(n2);
   Real* a = A.Store(); Real* b = B.Store(); Real* u = U.Store(); int i = n2;
   while (i--) { *a++ = *u++; *b++ = *u++; }
   FFT(A,B,A,B);
   int n21 = n2 + 1;
   X.ReDimension(n21); Y.ReDimension(n21);
   i = n2 - 1;
   a = A.Store(); b = B.Store();              // first els of A and B
   Real* an = a + i; Real* bn = b + i;        // last els of A and B
   Real* x = X.Store(); Real* y = Y.Store();  // first els of X and Y
   Real* xn = x + n2; Real* yn = y + n2;      // last els of X and Y

   *x++ = *a + *b; *y++ = 0.0;                // first complex element
   *xn-- = *a++ - *b++; *yn-- = 0.0;          // last complex element

   int j = -1; i = n2/2;
   while (i--)
   {
      Real c,s; cossin(j--,n,c,s);
      Real am = *a - *an; Real ap = *a++ + *an--;
      Real bm = *b - *bn; Real bp = *b++ + *bn--;
      Real samcbp = s * am + c * bp; Real sbpcam = s * bp - c * am;
      *x++  =  0.5 * ( ap + samcbp); *y++  =  0.5 * ( bm + sbpcam);
      *xn-- =  0.5 * ( ap - samcbp); *yn-- =  0.5 * (-bm + sbpcam);
   }
}

void RealFFTI(const ColumnVector& A, const ColumnVector& B, ColumnVector& U)
{
   // inverse of a Fourier transform of a real series
   Tracer trace("RealFFTI");
   const int n21 = A.Nrows();                     // length of arrays
   if (n21 != B.Nrows() || n21 == 0)
      Throw(ProgramException("Vector lengths unequal or zero", A, B));
   const int n2 = n21 - 1;  const int n = 2 * n2;  int i = n2 - 1;

   ColumnVector X(n2), Y(n2);
   Real* a = A.Store(); Real* b = B.Store();  // first els of A and B
   Real* an = a + n2;   Real* bn = b + n2;    // last els of A and B
   Real* x = X.Store(); Real* y = Y.Store();  // first els of X and Y
   Real* xn = x + i;    Real* yn = y + i;     // last els of X and Y

   Real hn = 0.5 / n2;
   *x++  = hn * (*a + *an);  *y++  = - hn * (*a - *an);
   a++; an--; b++; bn--;
   int j = -1;  i = n2/2;
   while (i--)
   {
      Real c,s; cossin(j--,n,c,s);
      Real am = *a - *an; Real ap = *a++ + *an--;
      Real bm = *b - *bn; Real bp = *b++ + *bn--;
      Real samcbp = s * am - c * bp; Real sbpcam = s * bp + c * am;
      *x++  =  hn * ( ap + samcbp); *y++  =  - hn * ( bm + sbpcam);
      *xn-- =  hn * ( ap - samcbp); *yn-- =  - hn * (-bm + sbpcam);
   }
   FFT(X,Y,X,Y);             // have done inverting elsewhere
   U.ReDimension(n); i = n2;
   x = X.Store(); y = Y.Store(); Real* u = U.Store();
   while (i--) { *u++ = *x++; *u++ = - *y++; }
}

