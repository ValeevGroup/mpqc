//$$ bandmat.cxx                     Band matrix definitions

// Copyright (C) 1991,2,3: R B Davies

#define WANT_MATH                    // include.h will get math fns

#include "include.h"

#include "newmat.h"
#include "newmatrc.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,4); ++ExeCount; }

#define REPORT {}

//#define REPORT1 { static ExeCounter ExeCount(__LINE__,4); ExeCount++; }

// REPORT1 constructors only - doesn't work in turbo and Borland C++

#define REPORT1 {}

//#define MONITOR(what,storage,store) \
//   { cout << what << " " << storage << " at " << (long)store << "\n"; }

#define MONITOR(what,store,storage) {}

BandMatrix::BandMatrix(const BaseMatrix& M)
{
   REPORT1 // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::BM);
   GetMatrix(gmx); CornerClear();
}

void BandMatrix::SetParameters(const GeneralMatrix* gmx)
{
   MatrixBandWidth bw = gmx->BandWidth();
   lower = bw.lower; upper = bw.upper;
}

void BandMatrix::ReDimension(int n, int lb, int ub)
{
   REPORT
   Tracer tr("BandMatrix::ReDimension");
   if (lb<0 || ub<0) Throw(ProgramException("Undefined bandwidth"));
   lower = (lb<=n) ? lb : n-1; upper = (ub<=n) ? ub : n-1;
   GeneralMatrix::ReDimension(n,n,n*(lower+1+upper)); CornerClear();
}

void BandMatrix::operator=(const BaseMatrix& X)
{
   REPORT // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::BM); CornerClear();
}

void BandMatrix::CornerClear() const
{
   // set unused parts of BandMatrix to zero
   REPORT
   int i = lower; Real* s = store; int bw = lower + 1 + upper;
   while (i)
      { int j = i--; Real* sj = s; s += bw; while (j--) *sj++ = 0.0; }
   i = upper; s = store + storage;
   while (i)
      { int j = i--; Real* sj = s; s -= bw; while (j--) *(--sj) = 0.0; }
}

MatrixBandWidth MatrixBandWidth::operator+(const MatrixBandWidth& bw) const
{
   int l = bw.lower; int u = bw.upper;
   l = (lower < 0 || l < 0) ? -1 : (lower > l) ? lower : l;
   u = (upper < 0 || u < 0) ? -1 : (upper > u) ? upper : u;
   return MatrixBandWidth(l,u);
}

MatrixBandWidth MatrixBandWidth::operator*(const MatrixBandWidth& bw) const
{
   int l = bw.lower; int u = bw.upper;
   l = (lower < 0 || l < 0) ? -1 : lower+l;
   u = (upper < 0 || u < 0) ? -1 : upper+u;
   return MatrixBandWidth(l,u);
}

UpperBandMatrix::UpperBandMatrix(const BaseMatrix& M)
{
   REPORT1 // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::UB);
   GetMatrix(gmx); CornerClear();
}

void UpperBandMatrix::operator=(const BaseMatrix& X)
{
   REPORT // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::UB); CornerClear();
}

LowerBandMatrix::LowerBandMatrix(const BaseMatrix& M)
{
   REPORT1 // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::LB);
   GetMatrix(gmx); CornerClear();
}

void LowerBandMatrix::operator=(const BaseMatrix& X)
{
   REPORT // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::LB); CornerClear();
}

BandLUMatrix::BandLUMatrix(const BaseMatrix& m)
{
   REPORT1
   Tracer tr("BandLUMatrix");
   GeneralMatrix* gm = ((BaseMatrix&)m).Evaluate(MatrixType::BM);
   GetMatrix(gm);
   m1 = ((BandMatrix*)gm)->lower; m2 = ((BandMatrix*)gm)->upper;
   if (nrows!=ncols) Throw(NotSquareException(*this));
   d = TRUE; sing = FALSE;
   indx = new int [nrows]; MatrixErrorNoSpace(indx);
   MONITOR_INT_NEW("Index (BndLUMat)",nrows,indx)
   storage2 = nrows * m1;
   store2 = new Real [storage2]; MatrixErrorNoSpace(store2);
   MONITOR_REAL_NEW("Make (BandLUMat)",storage2,store2)
   ludcmp();
}

BandLUMatrix::~BandLUMatrix()
{
   MONITOR_INT_DELETE("Index (BndLUMat)",nrows,indx)
   MONITOR_REAL_DELETE("Delete (BndLUMt)",storage2,store2)
#ifdef Version21
   delete [] indx; delete [] store2;
#else
   delete [nrows] indx; delete [storage2] store2;
#endif
}

MatrixType BandLUMatrix::Type() const { return MatrixType::BC; }


LogAndSign BandLUMatrix::LogDeterminant() const
{
   if (sing) return 0.0;
   Real* a = store; int w = m1+1+m2; LogAndSign sum; int i = nrows;
   while (i--) { sum *= *a; a += w; }
   if (!d) sum.ChangeSign(); return sum;
}

GeneralMatrix* BandMatrix::MakeSolver()
{
   REPORT
   GeneralMatrix* gm = new BandLUMatrix(*this);
   MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
}


void BandLUMatrix::ludcmp()
{
   REPORT
   Real* a = store;
   int i = m1; int j = m2; int k; int n = nrows; int w = m1 + 1 + m2;
   while (i)
   {
      Real* ai = a + i;
      k = ++j; while (k--) *a++ = *ai++;
      k = i--; while (k--) *a++ = 0.0;
   }

   a = store; int l = m1;
   for (k=0; k<n; k++)
   {
      Real x = *a; i = k; Real* aj = a;
      if (l < n) l++;
      for (j=k+1; j<l; j++)
         { aj += w; if (fabs(x) < fabs(*aj)) { x = *aj; i = j; } }
      indx[k] = i;
      if (x==0) { sing = TRUE; return; }
      if (i!=k)
      {
         d = !d; Real* ak = a; Real* ai = store + i * w; j = w;
         while (j--) { x = *ak; *ak++ = *ai; *ai++ = x; }
      }
      aj = a + w; Real* m = store2 + m1 * k;
      for (j=k+1; j<l; j++)
      {
         *m++ = x = *aj / *a; i = w; Real* ak = a;
	 while (--i) { Real* aj1 = aj++; *aj1 = *aj - x * *(++ak); }
         *aj++ = 0.0;
      }
      a += w;
   }
}

void BandLUMatrix::lubksb(Real* B, int mini)
{
   REPORT
   Tracer tr("BandLUMatrix::lubksb");
   if (sing) Throw(SingularException(*this));
   int n = nrows; int l = m1; int w = m1 + 1 + m2;

   for (int k=0; k<n; k++)
   {
      int i = indx[k];
      if (i!=k) { Real x=B[k]; B[k]=B[i]; B[i]=x; }
      if (l<n) l++;
      Real* m = store2 + k*m1; Real* b = B+k; Real* bi = b;
      for (i=k+1; i<l; i++)  *(++bi) -= *m++ * *b;
   }

   l = -m1;
   for (int i = n-1; i>=mini; i--)
   {
      Real* b = B + i; Real* bk = b; Real x = *bk;
      Real* a = store + w*i; Real y = *a;
      int k = l+m1; while (k--) x -=  *(++a) * *(++bk);
      *b = x / y;
      if (l < m2) l++;
   }
}

void BandLUMatrix::Solver(MatrixRowCol& mcout, const MatrixRowCol& mcin)
{
   REPORT
   Real* el = mcin.store; int i = mcin.skip;
   while (i--) *el++ = 0.0;
   el += mcin.storage; i = nrows - mcin.skip - mcin.storage;
   while (i--) *el++ = 0.0;
   lubksb(mcin.store, mcout.skip);
}

// Do we need check for entirely zero output?


void UpperBandMatrix::Solver(MatrixRowCol& mcout,
   const MatrixRowCol& mcin)
{
   REPORT
   Real* elx = mcin.store+mcout.skip; int i = mcin.skip-mcout.skip;
   while (i-- > 0) *elx++ = 0.0;
   int nr = mcin.skip+mcin.storage; elx = mcin.store+nr; Real* el = elx;
   int j = mcout.skip+mcout.storage-nr; i = nr-mcout.skip;
   while (j-- > 0) *elx++ = 0.0;

   Real* Ael = store + (upper+1)*(i-1)+1; j = 0;
   while (i-- > 0)
   {
      elx = el; Real sum = 0.0; int jx = j;
      while (jx--) sum += *(--Ael) * *(--elx);
      elx--; *elx = (*elx - sum) / *(--Ael);
      if (j<upper) Ael -= upper - (++j); else el--;
   }
}

void LowerBandMatrix::Solver(MatrixRowCol& mcout,
   const MatrixRowCol& mcin)
{
   REPORT
   Real* elx = mcin.store+mcout.skip; int i = mcin.skip-mcout.skip;
   while (i-- > 0) *elx++ = 0.0;
   int nc = mcin.skip; i = nc+mcin.storage; elx = mcin.store+i;
   int nr = mcout.skip+mcout.storage; int j = nr-i; i = nr-nc;
   while (j-- > 0) *elx++ = 0.0;

   Real* el = mcin.store+nc; Real* Ael = store + (lower+1)*nc + lower; j = 0;
   while (i-- > 0)
   {
      elx = el; Real sum = 0.0; int jx = j;
      while (jx--) sum += *Ael++ * *elx++;
      *elx = (*elx - sum) / *Ael++;
      if (j<lower) Ael += lower - (++j); else el++;
   }
}


LogAndSign BandMatrix::LogDeterminant() const
{
   REPORT
   BandLUMatrix C(*this); return C.LogDeterminant();
}

LogAndSign LowerBandMatrix::LogDeterminant() const
{
   REPORT
   int i = nrows; LogAndSign sum; Real* s = store + lower; int j = lower + 1;
   while (i--) { sum *= *s; s += j; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

LogAndSign UpperBandMatrix::LogDeterminant() const
{
   REPORT
   int i = nrows; LogAndSign sum; Real* s = store; int j = upper + 1;
   while (i--) { sum *= *s; s += j; }
   ((GeneralMatrix&)*this).tDelete(); return sum;
}

GeneralMatrix* SymmetricBandMatrix::MakeSolver()
{
   REPORT
   GeneralMatrix* gm = new BandLUMatrix(*this);
   MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
}

SymmetricBandMatrix::SymmetricBandMatrix(const BaseMatrix& M)
{
   REPORT1  // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::SB);
   GetMatrix(gmx);
}

GeneralMatrix* SymmetricBandMatrix::Transpose(TransposedMatrix*, MatrixType mt)
{ REPORT  return Evaluate(mt); }

LogAndSign SymmetricBandMatrix::LogDeterminant() const
{
   REPORT
   BandLUMatrix C(*this); return C.LogDeterminant();
}

void SymmetricBandMatrix::SetParameters(const GeneralMatrix* gmx)
{ lower = gmx->BandWidth().lower; }

void SymmetricBandMatrix::ReDimension(int n, int lb)
{
   REPORT
   Tracer tr("SymmetricBandMatrix::ReDimension");
   if (lb<0) Throw(ProgramException("Undefined bandwidth"));
   lower = (lb<=n) ? lb : n-1;
   GeneralMatrix::ReDimension(n,n,n*(lower+1));
}

void SymmetricBandMatrix::operator=(const BaseMatrix& X)
{
   REPORT // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::SB);
}

void SymmetricBandMatrix::CornerClear() const
{
   // set unused parts of BandMatrix to zero
   REPORT
   int i = lower; Real* s = store; int bw = lower + 1;
   while (i)
      { int j = i--; Real* sj = s; s += bw; while (j--) *sj++ = 0.0; }
}

MatrixBandWidth SymmetricBandMatrix::BandWidth() const
   { return MatrixBandWidth(lower,lower); }

inline Real square(Real x) { return x*x; }


Real SymmetricBandMatrix::SumSquare() const
{
   REPORT
   CornerClear();
   Real sum1=0.0; Real sum2=0.0; Real* s=store; int i=nrows; int l=lower;
   while (i--)
      { int j = l; while (j--) sum2 += square(*s++); sum1 += square(*s++); }
   ((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
}

Real SymmetricBandMatrix::SumAbsoluteValue() const
{
   REPORT
   CornerClear();
   Real sum1=0.0; Real sum2=0.0; Real* s=store; int i=nrows; int l=lower;
   while (i--)
      { int j = l; while (j--) sum2 += fabs(*s++); sum1 += fabs(*s++); }
   ((GeneralMatrix&)*this).tDelete(); return sum1 + 2.0 * sum2;
}

