//$$ newmat8.cxx         Advanced LU transform, scalar functions

// Copyright (C) 1991,2,3: R B Davies

#define WANT_MATH

#include "include.h"

#include "newmatap.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,8); ExeCount++; }

#define REPORT {}


/************************** LU transformation ****************************/

void CroutMatrix::ludcmp()
// LU decomposition - from numerical recipes in C
{
   REPORT
   Tracer trace("Crout(ludcmp)");
   int i,j;

   Real* vv=new Real [nrows]; MatrixErrorNoSpace(vv);
   MONITOR_REAL_NEW("Make  (CroutMat)",nrows,vv)
   Real* a;

   a=store;
   for (i=0;i<nrows;i++)
   {
      Real big=0.0;
      j=nrows; while (j--) { Real temp=fabs(*a++); if (temp > big) big=temp; }
      if (big == 0.0)
      {
         MONITOR_REAL_DELETE("Delete (CroutMt)",nrows,vv)
#ifdef Version21
         sing = TRUE; delete [] vv; return;
#else
         sing = TRUE; delete [nrows] vv; return;
#endif
      }
      vv[i]=1.0/big;
   }

   Real* aj=store;
   for (j=0;j<nrows;j++)
   {
      Real* ai=store;
      for (i=0;i<j;i++)
      {
         Real sum=*(ai+j); Real* aix=ai; Real* ajx=aj;
         int k=i; while (k--) { sum -= *aix++ * *ajx; ajx += nrows; }
         *ajx = sum; ai += nrows;
      }

      Real big=0.0; int imax;
      for (i=j;i<nrows;i++)
      {
         Real sum=*(ai+j); Real* aix=ai; Real* ajx=aj;
         int k=j; while (k--) { sum -= *aix++ * *ajx; ajx += nrows; }
         *aix = sum; ai += nrows;
         Real dum=vv[i]*fabs(sum); if (dum >= big) { big=dum; imax=i; }
      }

      if (j != imax)
      {
         Real* amax=store+imax*nrows; Real* ajx=store+j*nrows;
         int k=nrows; while(k--) { Real dum=*amax; *amax++=*ajx; *ajx++=dum; }
         d=!d; vv[imax]=vv[j];
      }

      indx[j]=imax; ai=aj+j*nrows;
      if (*ai == 0.0)
      {
         MONITOR_REAL_DELETE("Delete (CroutMt)",nrows,vv)
#ifdef Version21
         sing = TRUE; delete [] vv; return;
#else
         sing = TRUE; delete [nrows] vv; return;
#endif
      }
      Real dum=1.0/(*ai);
      i=nrows-j; while (--i) { ai += nrows; *ai *= dum; }

      aj++;
   }
   MONITOR_REAL_DELETE("Delete (CroutMt)",nrows,vv)
#ifdef Version21
   delete [] vv;
#else
   delete [nrows] vv;
#endif
}

void CroutMatrix::lubksb(Real* B, int mini)
{
   REPORT
   Tracer trace("Crout(lubksb)");
   if (sing) Throw(SingularException(*this));   
   int i,j; int ii=-1; Real* ai=store;

   for (i=0;i<nrows;i++)
   {
      int ip=indx[i]; Real sum=B[ip]; B[ip]=B[i];
      if (ii>=0)
      {
         Real* aij=ai+ii; Real* bj=B+ii; j=i-ii;
         while (j--) { sum -= *aij++ * *bj; bj++; }
      }
      else if (sum) ii=i;
      B[i]=sum; ai += nrows;
   }

   for (i=nrows-1;i>=mini;i--)
   {
      Real* bj=B+i; ai -= nrows; Real* ajx=ai+i; Real sum=*bj; bj++;
      j=nrows-i; while(--j) { sum -= *(++ajx) * *bj; bj++; }
      B[i]=sum/(*(ai+i));
   }
}


/****************************** scalar functions ****************************/

inline Real square(Real x) { return x*x; }

Real GeneralMatrix::SumSquare() const
{
   REPORT
   Real sum = 0.0; int i = storage; Real* s = store;
   while (i--) sum += square(*s++);
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real GeneralMatrix::SumAbsoluteValue() const
{
   REPORT
   Real sum = 0.0; int i = storage; Real* s = store;
   while (i--) sum += fabs(*s++);
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real GeneralMatrix::MaximumAbsoluteValue() const
{
   REPORT
   Real maxval = 0.0; int i = storage; Real* s = store;
   while (i--) { Real a = fabs(*s++); if (maxval < a) maxval = a; }
   ((GeneralMatrix&)*this).tDelete(); return maxval;
}

Real SymmetricMatrix::SumSquare() const
{
   REPORT
   Real sum1 = 0.0; Real sum2 = 0.0; Real* s = store; int nr = nrows;
   for (int i = 0; i<nr; i++)
   {
      int j = i;
      while (j--) sum2 += square(*s++);
      sum1 += square(*s++);
   }
   ((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
}

Real SymmetricMatrix::SumAbsoluteValue() const
{
   REPORT
   Real sum1 = 0.0; Real sum2 = 0.0; Real* s = store; int nr = nrows;
   for (int i = 0; i<nr; i++)
   {
      int j = i;
      while (j--) sum2 += fabs(*s++);
      sum1 += fabs(*s++);
   }
   ((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
}

Real BaseMatrix::SumSquare() const
{
   REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   Real s = gm->SumSquare(); return s;
}

Real BaseMatrix::SumAbsoluteValue() const
{
   REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   Real s = gm->SumAbsoluteValue(); return s;
}

Real BaseMatrix::MaximumAbsoluteValue() const
{
   REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   Real s = gm->MaximumAbsoluteValue(); return s;
}

Real Matrix::Trace() const
{
   REPORT
   Tracer trace("Trace");
   int i = nrows; int d = i+1;
   if (i != ncols) Throw(NotSquareException(*this));
   Real sum = 0.0; Real* s = store;
   while (i--) { sum += *s; s += d; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real DiagonalMatrix::Trace() const
{
   REPORT
   int i = nrows; Real sum = 0.0; Real* s = store;
   while (i--) sum += *s++;
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real SymmetricMatrix::Trace() const
{
   REPORT
   int i = nrows; Real sum = 0.0; Real* s = store; int j = 2;
   while (i--) { sum += *s; s += j++; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real LowerTriangularMatrix::Trace() const
{
   REPORT
   int i = nrows; Real sum = 0.0; Real* s = store; int j = 2;
   while (i--) { sum += *s; s += j++; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real UpperTriangularMatrix::Trace() const
{
   REPORT
   int i = nrows; Real sum = 0.0; Real* s = store;
   while (i) { sum += *s; s += i--; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real BandMatrix::Trace() const
{
   REPORT
   int i = nrows; int w = lower+upper+1;
   Real sum = 0.0; Real* s = store+lower;
   while (i--) { sum += *s; s += w; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real SymmetricBandMatrix::Trace() const
{
   REPORT
   int i = nrows; int w = lower+1;
   Real sum = 0.0; Real* s = store+lower;
   while (i--) { sum += *s; s += w; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

Real BaseMatrix::Trace() const
{
   REPORT
   GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate(MatrixType::Dg);
   Real sum = gm->Trace(); return sum;
}

void LogAndSign::operator*=(Real x)
{
   if (x > 0.0) { log_value += log(x); }
   else if (x < 0.0) { log_value += log(-x); sign = -sign; }
   else sign = 0;
}

Real LogAndSign::Value() const { return sign * exp(log_value); }

LogAndSign::LogAndSign(Real f)
{
   if (f == 0.0) { log_value = 0.0; sign = 0; return; }
   else if (f < 0.0) { sign = -1; f = -f; }
   else sign = 1;
   log_value = log(f);
}

LogAndSign DiagonalMatrix::LogDeterminant() const
{
   REPORT
   int i = nrows; LogAndSign sum; Real* s = store;
   while (i--) sum *= *s++;
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

LogAndSign LowerTriangularMatrix::LogDeterminant() const
{
   REPORT
   int i = nrows; LogAndSign sum; Real* s = store; int j = 2;
   while (i--) { sum *= *s; s += j++; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

LogAndSign UpperTriangularMatrix::LogDeterminant() const
{
   REPORT
   int i = nrows; LogAndSign sum; Real* s = store;
   while (i) { sum *= *s; s += i--; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

LogAndSign BaseMatrix::LogDeterminant() const
{
   REPORT GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   LogAndSign sum = gm->LogDeterminant(); return sum;
}

LogAndSign GeneralMatrix::LogDeterminant() const
{
   REPORT
   Tracer tr("Determinant");
   if (nrows != ncols) Throw(NotSquareException(*this));
   CroutMatrix C(*this); return C.LogDeterminant();
}

LogAndSign CroutMatrix::LogDeterminant() const
{
   REPORT
   if (sing) return 0.0;
   int i = nrows; int dd = i+1; LogAndSign sum; Real* s = store;
   while (i--) { sum *= *s; s += dd; }
   if (!d) sum.ChangeSign(); return sum;

}

LinearEquationSolver::LinearEquationSolver(const BaseMatrix& bm)
: gm( ( ((BaseMatrix&)bm).Evaluate() )->MakeSolver() )
{
   if (gm==&bm) { REPORT  gm = gm->Image(); } 
   // want a copy if  *gm is actually bm
   else { REPORT  gm->Protect(); }
}

