//$$ newmat6.cxx            Operators, element access, submatrices

// Copyright (C) 1991,2,3: R B Davies

#define WANT_STREAM
#include "include.h"

#include "newmat.h"
#include "newmatrc.h"


//#define REPORT { static ExeCounter ExeCount(__LINE__,6); ExeCount++; }

#define REPORT {}

/*************************** general utilities *************************/

static int tristore(int n)                      // els in triangular matrix
{ return (n*(n+1))/2; }


/****************************** operators *******************************/

Real& Matrix::operator()(int m, int n)
{
  if (m<=0 || m>nrows || n<=0 || n>ncols){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
  return store[(m-1)*ncols+n-1];
}

Real& SymmetricMatrix::operator()(int m, int n)
{
  if (m<=0 || n<=0 || m>nrows || n>ncols){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
   if (m>=n) return store[tristore(m-1)+n-1];
   else return store[tristore(n-1)+m-1];
}

Real& UpperTriangularMatrix::operator()(int m, int n)
{
  if (m<=0 || n<m || n>ncols){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
  return store[(m-1)*ncols+n-1-tristore(m-1)];
}

Real& LowerTriangularMatrix::operator()(int m, int n)
{
  if (n<=0 || m<n || m>nrows){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
   return store[tristore(m-1)+n-1];
}

Real& DiagonalMatrix::operator()(int m, int n)
{
  if (n<=0 || m!=n || m>nrows || n>ncols){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
  return store[n-1];
}

Real& DiagonalMatrix::operator()(int m)
{
   if (m<=0 || m>nrows) {
     cerr << "Subscript out of bounds\n";
     Throw(IndexException(m,*this));
   }
   return store[m-1];
}

Real& ColumnVector::operator()(int m)
{
   if (m<=0 || m> nrows) {
     cerr << "Subscript out of bounds\n";
     Throw(IndexException(m,*this));
   }
   return store[m-1];
}

Real& RowVector::operator()(int n)
{
   if (n<=0 || n> ncols) {
     cerr << "Subscript out of bounds\n";
     Throw(IndexException(n,*this));
   }
   return store[n-1];
}

Real& BandMatrix::operator()(int m, int n)
{
  int w = upper+lower+1; int i = lower+n-m;
  if (m<=0 || m>nrows || n<=0 || n>ncols || i<0 || i>=w){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
  return store[w*(m-1)+i];
}

Real& UpperBandMatrix::operator()(int m, int n)
{
  int w = upper+1; int i = n-m;
  if (m<=0 || m>nrows || n<=0 || n>ncols || i<0 || i>=w){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
  return store[w*(m-1)+i];
}

Real& LowerBandMatrix::operator()(int m, int n)
{
  int w = lower+1; int i = lower+n-m;
  if (m<=0 || m>nrows || n<=0 || n>ncols || i<0 || i>=w){
    cerr << "Subscript out of bounds\n";
    Throw(IndexException(m,n,*this));
  }
  return store[w*(m-1)+i];
}

Real& SymmetricBandMatrix::operator()(int m, int n)
{
   int w = lower+1;
   if (m>=n)
   {
     int i = lower+n-m;
     if ( m>nrows || n<=0 || i<0 ){
       cerr << "Subscript out of bounds\n";
       Throw(IndexException(m,n,*this));
     }
     return store[w*(m-1)+i];
   }
   else
   {
     int i = lower+m-n;
     if ( n>nrows || m<=0 || i<0 ){
       cerr << "Subscript out of bounds\n";
       Throw(IndexException(m,n,*this));
     }
      return store[w*(n-1)+i];
   }
}

#ifndef __ZTC__

Real Matrix::operator()(int m, int n) const
{
   if (m<=0 || m>nrows || n<=0 || n>ncols)
      Throw(IndexException(m,n,*this));
   return store[(m-1)*ncols+n-1];
}

Real SymmetricMatrix::operator()(int m, int n) const
{
   if (m<=0 || n<=0 || m>nrows || n>ncols)
      Throw(IndexException(m,n,*this));
   if (m>=n) return store[tristore(m-1)+n-1];
   else return store[tristore(n-1)+m-1];
}

Real UpperTriangularMatrix::operator()(int m, int n) const
{
   if (m<=0 || n<m || n>ncols)
      Throw(IndexException(m,n,*this));
   return store[(m-1)*ncols+n-1-tristore(m-1)];
}

Real LowerTriangularMatrix::operator()(int m, int n) const
{
   if (n<=0 || m<n || m>nrows)
      Throw(IndexException(m,n,*this));
   return store[tristore(m-1)+n-1];
}

Real DiagonalMatrix::operator()(int m, int n) const
{
   if (n<=0 || m!=n || m>nrows || n>ncols)
      Throw(IndexException(m,n,*this));
   return store[n-1];
}

Real DiagonalMatrix::operator()(int m) const
{
   if (m<=0 || m>nrows) Throw(IndexException(m,*this));
   return store[m-1];
}

Real ColumnVector::operator()(int m) const
{
   if (m<=0 || m> nrows) Throw(IndexException(m,*this));
   return store[m-1];
}

Real RowVector::operator()(int n) const
{
   if (n<=0 || n> ncols) Throw(IndexException(n,*this));
   return store[n-1];
}

Real BandMatrix::operator()(int m, int n) const
{
   int w = upper+lower+1; int i = lower+n-m;
   if (m<=0 || m>nrows || n<=0 || n>ncols || i<0 || i>=w)
      Throw(IndexException(m,n,*this));
   return store[w*(m-1)+i];
}

Real SymmetricBandMatrix::operator()(int m, int n) const
{
   int w = lower+1;
   if (m>=n)
   {
      int i = lower+n-m;
      if ( m>nrows || n<=0 || i<0 )
         Throw(IndexException(m,n,*this));
      return store[w*(m-1)+i];
   }
   else
   {
      int i = lower+m-n;
      if ( n>nrows || m<=0 || i<0 )
         Throw(IndexException(m,n,*this));
      return store[w*(n-1)+i];
   }
}

#endif

Real BaseMatrix::AsScalar() const
{
   REPORT
   GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   if (gm->nrows!=1 || gm->ncols!=1)
   {
      Try
      {
         Tracer tr("AsScalar");
	 Throw(ProgramException("Cannot convert to scalar", *gm));
      }
      CatchAll { gm->tDelete(); Throw(); }
   }
   Real x = *(gm->store); gm->tDelete(); return x;
}

#ifdef TEMPS_DESTROYED_QUICKLY

AddedMatrix& BaseMatrix::operator+(const BaseMatrix& bm) const
{
   REPORT
   AddedMatrix* x = new AddedMatrix(this, &bm);
   MatrixErrorNoSpace(x);
   return *x;
}

MultipliedMatrix& BaseMatrix::operator*(const BaseMatrix& bm) const
{
   REPORT
   MultipliedMatrix* x = new MultipliedMatrix(this, &bm);
   MatrixErrorNoSpace(x);
   return *x;
}

//SolvedMatrix& InvertedMatrix::operator*(const BaseMatrix& bmx) const
SolvedMatrix& InvertedMatrix::operator*(const BaseMatrix& bmx)
{
   REPORT
   SolvedMatrix* x;
   Try { x = new SolvedMatrix(bm, &bmx); MatrixErrorNoSpace(x); }
   CatchAll { delete this; Throw(); }
   delete this;                // since we are using bm rather than this
   return *x;
}

SubtractedMatrix& BaseMatrix::operator-(const BaseMatrix& bm) const
{
   REPORT
   SubtractedMatrix* x = new SubtractedMatrix(this, &bm);
   MatrixErrorNoSpace(x);
   return *x;
}

ShiftedMatrix& BaseMatrix::operator+(Real f) const
{
   REPORT
   ShiftedMatrix* x = new ShiftedMatrix(this, f);
   MatrixErrorNoSpace(x);
   return *x;
}

ScaledMatrix& BaseMatrix::operator*(Real f) const
{
   REPORT
   ScaledMatrix* x = new ScaledMatrix(this, f);
   MatrixErrorNoSpace(x);
   return *x;
}

ScaledMatrix& BaseMatrix::operator/(Real f) const
{
   REPORT
   ScaledMatrix* x = new ScaledMatrix(this, 1.0/f);
   MatrixErrorNoSpace(x);
   return *x;
}

ShiftedMatrix& BaseMatrix::operator-(Real f) const
{
   REPORT
   ShiftedMatrix* x = new ShiftedMatrix(this, -f);
   MatrixErrorNoSpace(x);
   return *x;
}

TransposedMatrix& BaseMatrix::t() const
{
   REPORT
   TransposedMatrix* x = new TransposedMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

NegatedMatrix& BaseMatrix::operator-() const
{
   REPORT
   NegatedMatrix* x = new NegatedMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

InvertedMatrix& BaseMatrix::i() const
{
   REPORT
   InvertedMatrix* x = new InvertedMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

ConstMatrix& GeneralMatrix::c() const
{
   if (tag != -1)
      Throw(ProgramException(".c() applied to temporary matrix"));
   REPORT
   ConstMatrix* x = new ConstMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

RowedMatrix& BaseMatrix::AsRow() const
{
   REPORT
   RowedMatrix* x = new RowedMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

ColedMatrix& BaseMatrix::AsColumn() const
{
   REPORT
   ColedMatrix* x = new ColedMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

DiagedMatrix& BaseMatrix::AsDiagonal() const
{
   REPORT
   DiagedMatrix* x = new DiagedMatrix(this);
   MatrixErrorNoSpace(x);
   return *x;
}

MatedMatrix& BaseMatrix::AsMatrix(int nrx, int ncx) const
{
   REPORT
   MatedMatrix* x = new MatedMatrix(this,nrx,ncx);
   MatrixErrorNoSpace(x);
   return *x;
}

#else

AddedMatrix BaseMatrix::operator+(const BaseMatrix& bm) const
{ REPORT return AddedMatrix(this, &bm); }

MultipliedMatrix BaseMatrix::operator*(const BaseMatrix& bm) const
{ REPORT return MultipliedMatrix(this, &bm); }

SolvedMatrix InvertedMatrix::operator*(const BaseMatrix& bmx) const
{ REPORT return SolvedMatrix(bm, &bmx); }

SubtractedMatrix BaseMatrix::operator-(const BaseMatrix& bm) const
{ REPORT return SubtractedMatrix(this, &bm); }

ShiftedMatrix BaseMatrix::operator+(Real f) const
{ REPORT return ShiftedMatrix(this, f); }

ScaledMatrix BaseMatrix::operator*(Real f) const
{ REPORT return ScaledMatrix(this, f); }

ScaledMatrix BaseMatrix::operator/(Real f) const
{ REPORT return ScaledMatrix(this, 1.0/f); }

ShiftedMatrix BaseMatrix::operator-(Real f) const
{ REPORT return ShiftedMatrix(this, -f); }

TransposedMatrix BaseMatrix::t() const
{ REPORT return TransposedMatrix(this); }

NegatedMatrix BaseMatrix::operator-() const
{ REPORT return NegatedMatrix(this); }

InvertedMatrix BaseMatrix::i() const
{ REPORT return InvertedMatrix(this); }

ConstMatrix GeneralMatrix::c() const
{
   if (tag != -1)
      Throw(ProgramException(".c() applied to temporary matrix"));
   REPORT return ConstMatrix(this);
}

RowedMatrix BaseMatrix::AsRow() const
{ REPORT return RowedMatrix(this); }

ColedMatrix BaseMatrix::AsColumn() const
{ REPORT return ColedMatrix(this); }

DiagedMatrix BaseMatrix::AsDiagonal() const
{ REPORT return DiagedMatrix(this); }

MatedMatrix BaseMatrix::AsMatrix(int nrx, int ncx) const
{ REPORT return MatedMatrix(this,nrx,ncx); }

#endif

void GeneralMatrix::operator=(Real f)
{ REPORT int i=storage; Real* s=store; while (i--) { *s++ = f; } }

void Matrix::operator=(const BaseMatrix& X)
{
   REPORT //CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::Rt);
} 

void RowVector::operator=(const BaseMatrix& X)
{
   REPORT  // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::RV);
   if (nrows!=1)
   {
      Tracer tr("RowVector(=)");
      Throw(VectorException(*this));
   }
}

void ColumnVector::operator=(const BaseMatrix& X)
{
   REPORT //CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::CV);
   if (ncols!=1)
   {
      Tracer tr("ColumnVector(=)");
      Throw(VectorException(*this));
   }
}

void SymmetricMatrix::operator=(const BaseMatrix& X)
{
   REPORT // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::Sm);
}
 
void UpperTriangularMatrix::operator=(const BaseMatrix& X)
{
   REPORT //CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::UT);
}

void LowerTriangularMatrix::operator=(const BaseMatrix& X)
{
   REPORT //CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::LT);
}

void DiagonalMatrix::operator=(const BaseMatrix& X)
{
   REPORT // CheckConversion(X);
   MatrixConversionCheck mcc;
   Eq(X,MatrixType::Dg);
}

void GeneralMatrix::operator<<(const Real* r)
{
   REPORT
   int i = storage; Real* s=store;
   while(i--) *s++ = *r++;
}


/************************* element access *********************************/

Real& Matrix::element(int m, int n)
{
   if (m<0 || m>= nrows || n<0 || n>= ncols)
      Throw(IndexException(m,n,*this,TRUE));
   return store[m*ncols+n];
}

const Real& Matrix::element(int m, int n) const
{
   if (m<0 || m>= nrows || n<0 || n>= ncols)
      Throw(IndexException(m,n,*this,TRUE));
   return store[m*ncols+n];
}

Real& SymmetricMatrix::element(int m, int n)
{
   if (m<0 || n<0 || m >= nrows || n>=ncols)
      Throw(IndexException(m,n,*this,TRUE));
   if (m>=n) return store[tristore(m)+n];
   else return store[tristore(n)+m];
}

const Real& SymmetricMatrix::element(int m, int n) const
{
   if (m<0 || n<0 || m >= nrows || n>=ncols)
      Throw(IndexException(m,n,*this,TRUE));
   if (m>=n) return store[tristore(m)+n];
   else return store[tristore(n)+m];
}

Real& UpperTriangularMatrix::element(int m, int n)
{
   if (m<0 || n<m || n>=ncols)
      Throw(IndexException(m,n,*this,TRUE));
   return store[m*ncols+n-tristore(m)];
}

Real& LowerTriangularMatrix::element(int m, int n)
{
   if (n<0 || m<n || m>=nrows)
      Throw(IndexException(m,n,*this,TRUE));
   return store[tristore(m)+n];
}

Real& DiagonalMatrix::element(int m, int n)
{
   if (n<0 || m!=n || m>=nrows || n>=ncols)
      Throw(IndexException(m,n,*this,TRUE));
   return store[n];
}

Real& DiagonalMatrix::element(int m)
{
   if (m<0 || m>=nrows) Throw(IndexException(m,*this,TRUE));
   return store[m];
}

const Real& DiagonalMatrix::element(int m, int n) const
{
   if (n<0 || m!=n || m>=nrows || n>=ncols)
      Throw(IndexException(m,n,*this,TRUE));
   return store[n];
}

const Real& DiagonalMatrix::element(int m) const
{
   if (m<0 || m>=nrows) Throw(IndexException(m,*this,TRUE));
   return store[m];
}

Real& ColumnVector::element(int m)
{
   if (m<0 || m>= nrows) Throw(IndexException(m,*this,TRUE));
   return store[m];
}

const Real& ColumnVector::element(int m) const
{
   if (m<0 || m>= nrows) Throw(IndexException(m,*this,TRUE));
   return store[m];
}

Real& RowVector::element(int n)
{
   if (n<0 || n>= ncols)  Throw(IndexException(n,*this,TRUE));
   return store[n];
}

const Real& RowVector::element(int n) const
{
   if (n<0 || n>= ncols)  Throw(IndexException(n,*this,TRUE));
   return store[n];
}

Real& BandMatrix::element(int m, int n)
{
   int w = upper+lower+1; int i = lower+n-m;
   if (m<0 || m>= nrows || n<0 || n>= ncols || i<0 || i>=w)
      Throw(IndexException(m,n,*this,TRUE));
   return store[w*m+i];
}

Real& UpperBandMatrix::element(int m, int n)
{
   int w = upper+1; int i = n-m;
   if (m<0 || m>= nrows || n<0 || n>= ncols || i<0 || i>=w)
      Throw(IndexException(m,n,*this,TRUE));
   return store[w*m+i];
}

Real& LowerBandMatrix::element(int m, int n)
{
   int w = lower+1; int i = lower+n-m;
   if (m<0 || m>= nrows || n<0 || n>= ncols || i<0 || i>=w)
      Throw(IndexException(m,n,*this,TRUE));
   return store[w*m+i];
}

Real& SymmetricBandMatrix::element(int m, int n)
{
   int w = lower+1;
   if (m>=n)
   {
      int i = lower+n-m;
      if ( m>=nrows || n<0 || i<0 )
         Throw(IndexException(m,n,*this,TRUE));
      return store[w*m+i];
   }
   else
   {
      int i = lower+m-n;
      if ( n>=nrows || m<0 || i<0 )
         Throw(IndexException(m,n,*this,TRUE));
      return store[w*n+i];
   }
}

