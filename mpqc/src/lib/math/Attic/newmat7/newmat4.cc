//$$ newmat4.cxx       Constructors, ReDimension, basic utilities

// Copyright (C) 1991,2,3: R B Davies

#include "include.h"

#include "newmat.h"
#include "newmatrc.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,4); ExeCount++; }

#define REPORT {}

//#define REPORT1 { static ExeCounter ExeCount(__LINE__,4); ExeCount++; }

// REPORT1 constructors only - doesn't work in turbo and Borland C++

#define REPORT1 {}


/*************************** general utilities *************************/

static int tristore(int n)                      // els in triangular matrix
{ return (n*(n+1))/2; }


/****************************** constructors ***************************/

GeneralMatrix::GeneralMatrix()
{ store=0; storage=0; nrows=0; ncols=0; tag=-1; }

GeneralMatrix::GeneralMatrix(ArrayLengthSpecifier s)
{
   REPORT1
   storage=s.Value(); tag=-1;
   if (storage)
   {
      store = new Real [storage]; MatrixErrorNoSpace(store);
      MONITOR_REAL_NEW("Make (GenMatrix)",storage,store)
   }
   else store = 0;
}

Matrix::Matrix(int m, int n) : GeneralMatrix(m*n)
{ REPORT1 nrows=m; ncols=n; }

SymmetricMatrix::SymmetricMatrix(ArrayLengthSpecifier n)
   : GeneralMatrix(tristore(n.Value()))
{ REPORT1 nrows=n.Value(); ncols=n.Value(); }

UpperTriangularMatrix::UpperTriangularMatrix(ArrayLengthSpecifier n)
   : GeneralMatrix(tristore(n.Value()))
{ REPORT1 nrows=n.Value(); ncols=n.Value(); }

LowerTriangularMatrix::LowerTriangularMatrix(ArrayLengthSpecifier n)
   : GeneralMatrix(tristore(n.Value()))
{ REPORT1 nrows=n.Value(); ncols=n.Value(); }

DiagonalMatrix::DiagonalMatrix(ArrayLengthSpecifier m) : GeneralMatrix(m)
{ REPORT1 nrows=m.Value(); ncols=m.Value(); }

Matrix::Matrix(const BaseMatrix& M)
{
   REPORT1 // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Rt);
   GetMatrix(gmx);
}

RowVector::RowVector(const BaseMatrix& M) : Matrix(M)
{
   if (nrows!=1)
   {
      Tracer tr("RowVector");
      Throw(VectorException(*this));
   }
}

ColumnVector::ColumnVector(const BaseMatrix& M) : Matrix(M)
{
   if (ncols!=1)
   {
      Tracer tr("ColumnVector");
      Throw(VectorException(*this));
   }
}

SymmetricMatrix::SymmetricMatrix(const BaseMatrix& M)
{
   REPORT1  // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Sm);
   GetMatrix(gmx);
}

UpperTriangularMatrix::UpperTriangularMatrix(const BaseMatrix& M)
{
   REPORT1 // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::UT);
   GetMatrix(gmx);
}

LowerTriangularMatrix::LowerTriangularMatrix(const BaseMatrix& M)
{
   REPORT1 // CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::LT);
   GetMatrix(gmx);
}

DiagonalMatrix::DiagonalMatrix(const BaseMatrix& M)
{
   REPORT1 //CheckConversion(M);
   MatrixConversionCheck mcc;
   GeneralMatrix* gmx=((BaseMatrix&)M).Evaluate(MatrixType::Dg);
   GetMatrix(gmx);
}

GeneralMatrix::~GeneralMatrix()
{
   if (store)
   {
      MONITOR_REAL_DELETE("Free (GenMatrix)",storage,store)
#ifdef Version21
      delete [] store;
#else
      delete [storage] store;
#endif
   }
}

CroutMatrix::CroutMatrix(const BaseMatrix& m)
{
   REPORT1
   Tracer tr("CroutMatrix");
   GeneralMatrix* gm = ((BaseMatrix&)m).Evaluate(MatrixType::Rt);
   GetMatrix(gm);
   if (nrows!=ncols) Throw(NotSquareException(*this));
   d=TRUE; sing=FALSE;
   indx=new int [nrows]; MatrixErrorNoSpace(indx);
   MONITOR_INT_NEW("Index (CroutMat)",nrows,indx)
   ludcmp();
}

CroutMatrix::~CroutMatrix()
{
   MONITOR_INT_DELETE("Index (CroutMat)",nrows,indx)
#ifdef Version21
   delete [] indx;
#else
   delete [nrows] indx;
#endif
}

//ReturnMatrixX::ReturnMatrixX(GeneralMatrix& gmx)
//{
//   REPORT1
//   gm = gmx.Image(); gm->ReleaseAndDelete();
//}

#ifndef TEMPS_DESTROYED_QUICKLY

GeneralMatrix::operator ReturnMatrixX() const
{
   REPORT
   GeneralMatrix* gm = Image(); gm->ReleaseAndDelete(); 
   return ReturnMatrixX(gm);
}

#else

GeneralMatrix::operator ReturnMatrixX&() const
{
   REPORT
   GeneralMatrix* gm = Image(); gm->ReleaseAndDelete();
   ReturnMatrixX* x = new ReturnMatrixX(gm);
   MatrixErrorNoSpace(x); return *x;
}

#endif

/**************************** ReDimension matrices ***************************/

void GeneralMatrix::ReDimension(int nr, int nc, int s)
{
   REPORT 
   if (store)
   {
      MONITOR_REAL_DELETE("Free (ReDimensi)",storage,store)
#ifdef Version21
      delete [] store;
#else
      delete [storage] store;
#endif
   }
   storage=s; nrows=nr; ncols=nc; tag=-1;
   if (s)
   {
      store = new Real [storage]; MatrixErrorNoSpace(store);
      MONITOR_REAL_NEW("Make (ReDimensi)",storage,store)
   }
   else store = 0;
}

void Matrix::ReDimension(int nr, int nc)
{ REPORT GeneralMatrix::ReDimension(nr,nc,nr*nc); }

void SymmetricMatrix::ReDimension(int nr)
{ REPORT GeneralMatrix::ReDimension(nr,nr,tristore(nr)); }

void UpperTriangularMatrix::ReDimension(int nr)
{ REPORT GeneralMatrix::ReDimension(nr,nr,tristore(nr)); }

void LowerTriangularMatrix::ReDimension(int nr)
{ REPORT GeneralMatrix::ReDimension(nr,nr,tristore(nr)); }

void DiagonalMatrix::ReDimension(int nr)
{ REPORT GeneralMatrix::ReDimension(nr,nr,nr); }

void RowVector::ReDimension(int nc)
{ REPORT GeneralMatrix::ReDimension(1,nc,nc); }

void ColumnVector::ReDimension(int nr)
{ REPORT GeneralMatrix::ReDimension(nr,1,nr); }


/********************* manipulate types, storage **************************/

int GeneralMatrix::search(const BaseMatrix* s) const
{ REPORT return (s==this) ? 1 : 0; }

int MultipliedMatrix::search(const BaseMatrix* s) const
{ REPORT return bm1->search(s) + bm2->search(s); }

int ShiftedMatrix::search(const BaseMatrix* s) const
{ REPORT return bm->search(s); }

int NegatedMatrix::search(const BaseMatrix* s) const
{ REPORT return bm->search(s); }

int ConstMatrix::search(const BaseMatrix* s) const
{ REPORT return (s==cgm) ? 1 : 0; }

int ReturnMatrixX::search(const BaseMatrix* s) const
{ REPORT return (s==gm) ? 1 : 0; }

MatrixType Matrix::Type() const { return MatrixType::Rt; }
MatrixType SymmetricMatrix::Type() const { return MatrixType::Sm; }
MatrixType UpperTriangularMatrix::Type() const { return MatrixType::UT; }
MatrixType LowerTriangularMatrix::Type() const { return MatrixType::LT; }
MatrixType DiagonalMatrix::Type() const { return MatrixType::Dg; }
MatrixType RowVector::Type() const { return MatrixType::RV; }
MatrixType ColumnVector::Type() const { return MatrixType::CV; }
MatrixType CroutMatrix::Type() const { return MatrixType::Ct; }
MatrixType BandMatrix::Type() const { return MatrixType::BM; }
MatrixType UpperBandMatrix::Type() const { return MatrixType::UB; }
MatrixType LowerBandMatrix::Type() const { return MatrixType::LB; }
MatrixType SymmetricBandMatrix::Type() const { return MatrixType::SB; }

MatrixBandWidth BaseMatrix::BandWidth() const { return -1; }
MatrixBandWidth DiagonalMatrix::BandWidth() const { return 0; }

MatrixBandWidth BandMatrix::BandWidth() const
   { return MatrixBandWidth(lower,upper); }

MatrixBandWidth AddedMatrix::BandWidth() const
{ return gm1->BandWidth() + gm2->BandWidth(); }

MatrixBandWidth MultipliedMatrix::BandWidth() const
{ return gm1->BandWidth() * gm2->BandWidth(); }

MatrixBandWidth SolvedMatrix::BandWidth() const { return -1; }
MatrixBandWidth ScaledMatrix::BandWidth() const { return gm->BandWidth(); }
MatrixBandWidth NegatedMatrix::BandWidth() const { return gm->BandWidth(); }

MatrixBandWidth TransposedMatrix::BandWidth() const
{ return gm->BandWidth().t(); }

MatrixBandWidth InvertedMatrix::BandWidth() const { return -1; }
MatrixBandWidth RowedMatrix::BandWidth() const { return -1; }
MatrixBandWidth ColedMatrix::BandWidth() const { return -1; }
MatrixBandWidth DiagedMatrix::BandWidth() const { return 0; }
MatrixBandWidth MatedMatrix::BandWidth() const { return -1; }
MatrixBandWidth ConstMatrix::BandWidth() const { return cgm->BandWidth(); }
MatrixBandWidth ReturnMatrixX::BandWidth() const { return gm->BandWidth(); }

MatrixBandWidth GetSubMatrix::BandWidth() const
{

   if (row_skip==col_skip && row_number==col_number) return gm->BandWidth();
   else return MatrixBandWidth(-1);
}

/************************ the memory managment tools **********************/

//  Rules regarding tDelete, reuse, GetStore
//    All matrices processed during expression evaluation must be subject
//    to exactly one of reuse(), tDelete(), GetStore() or BorrowStore().
//    If reuse returns TRUE the matrix must be reused.
//    GetMatrix(gm) always calls gm->GetStore()
//    gm->Evaluate obeys rules
//    bm->Evaluate obeys rules for matrices in bm structure

void GeneralMatrix::tDelete()
{
   if (tag<0)
   {
      if (tag<-1) { REPORT store=0; delete this; return; }  // borrowed
      else { REPORT return; }
   }
   if (tag==1)
   {
      REPORT  MONITOR_REAL_DELETE("Free   (tDelete)",storage,store)
#ifdef Version21
      if (store) delete [] store;
#else
      if (store) delete [storage] store;
#endif
      store=0; tag=-1; return;
   }
   if (tag==0) { REPORT delete this; return; }
   REPORT tag--; return;
}

static void BlockCopy(int n, Real* from, Real* to)
{
   REPORT
   int i = (n >> 3);
   while (i--)
   {
      *to++ = *from++; *to++ = *from++; *to++ = *from++; *to++ = *from++;
      *to++ = *from++; *to++ = *from++; *to++ = *from++; *to++ = *from++;
   }
   i = n & 7; while (i--) *to++ = *from++;
}

Boolean GeneralMatrix::reuse()
{
   if (tag<-1)
   {
      REPORT
      Real* s = new Real [storage]; MatrixErrorNoSpace(s);
      MONITOR_REAL_NEW("Make     (reuse)",storage,s)
      BlockCopy(storage, store, s); store=s; tag=0; return TRUE;
   }
   if (tag<0) { REPORT return FALSE; }
   if (tag<=1)  { REPORT return TRUE; }
   REPORT tag--; return FALSE;
}

Real* GeneralMatrix::GetStore()
{
   if (tag<0 || tag>1)
   {
      Real* s;
      if (storage)
      {
         s = new Real [storage]; MatrixErrorNoSpace(s);
         MONITOR_REAL_NEW("Make  (GetStore)",storage,s)
         BlockCopy(storage, store, s);
      }
      else s = 0;
      if (tag>1) { REPORT tag--; }
      else if (tag < -1) { REPORT store=0; delete this; } // borrowed store
      else { REPORT }
      return s;
   }
   Real* s=store; store=0;
   if (tag==0) { REPORT delete this; }
   else { REPORT tag=-1; }
   return s;
}

/*
#ifndef __ZTC__
void GeneralMatrix::GetMatrixC(const GeneralMatrix* gmx)
{
   REPORT tag=-1;
   nrows=gmx->nrows; ncols=gmx->ncols; storage=gmx->storage;
   SetParameters(gmx); 
   store = new Real [storage]; MatrixErrorNoSpace(store);
   MONITOR_REAL_NEW("Make (GetMatrix)",storage,store)
   BlockCopy(storage, gmx->store, store);
}
#endif
*/

void GeneralMatrix::GetMatrix(const GeneralMatrix* gmx)
{
   REPORT  tag=-1; nrows=gmx->Nrows(); ncols=gmx->Ncols();
   storage=gmx->storage; SetParameters(gmx);
   store=((GeneralMatrix*)gmx)->GetStore();
}

GeneralMatrix* GeneralMatrix::BorrowStore(GeneralMatrix* gmx, MatrixType mt)
// Copy storage of *this to storage of *gmx. Then convert to type mt.
// If mt == 0 just let *gmx point to storage of *this if tag==-1.
{
   if (!mt)
   {
      if (tag == -1) { REPORT gmx->tag = -2; gmx->store = store; }
      else { REPORT gmx->tag = 0; gmx->store = GetStore(); }
   }
   else if (Compare(gmx->Type(),mt))
   { REPORT  gmx->tag = 0; gmx->store = GetStore(); }
   else
   {
      REPORT gmx->tag = -2; gmx->store = store;
      gmx = gmx->Evaluate(mt); gmx->tag = 0; tDelete();
   }

   return gmx;
}

void GeneralMatrix::Eq(const BaseMatrix& X, MatrixType mt)
// Count number of references to this in X.
// If zero delete storage in X;
// otherwise tag X to show when storage can be deleted
// evaluate X and copy to current object
{
   int counter=X.search(this);
   if (counter==0)
   {
      REPORT
      if (store)
      {
         MONITOR_REAL_DELETE("Free (operator=)",storage,store)
#ifdef Version21
         REPORT delete [] store; storage=0;
#else
         REPORT delete [storage] store; storage=0;
#endif
      }
   }
   else { REPORT Release(counter); }
   GeneralMatrix* gmx = ((BaseMatrix&)X).Evaluate(mt);
   if (gmx!=this) { REPORT GetMatrix(gmx); }
   else { REPORT }
   Protect();
}

void GeneralMatrix::Inject(const GeneralMatrix& X)
// copy stored values of X; otherwise leave els of *this unchanged
{
   REPORT
   Tracer tr("Inject");
   if (nrows != X.nrows || ncols != X.ncols)
      Throw(IncompatibleDimensionsException());
   MatrixRow mr((GeneralMatrix*)&X, LoadOnEntry);
   MatrixRow mrx(this, LoadOnEntry+StoreOnExit+DirectPart);
   int i=nrows;
   while (i--) { mrx.Inject(mr); mrx.Next(); mr.Next(); }
}  

/*************** checking for data loss during conversion *******************/

//void GeneralMatrix::CheckConversion(const BaseMatrix& M)
//{
//   if (!(this->Type() >= M.Type()))
//      Throw(ProgramException("Illegal Conversion"));
//}

Boolean MatrixConversionCheck::DoCheck = FALSE;

void MatrixConversionCheck::DataLoss()
   { if (DoCheck) Throw(ProgramException("Illegal Conversion")); }

Boolean Compare(const MatrixType& source, MatrixType& destination)
{
   if (!destination) { destination=source; return TRUE; }
   if (destination==source) return TRUE;
   if (MatrixConversionCheck::IsOn() && !(destination>=source))
      Throw(ProgramException("Illegal Conversion"));
   return FALSE;
}

// Added my Mike Colvin
Boolean Rectangular(MatrixType a, MatrixType b, MatrixType c)
{ return ((a.attribute | b.attribute | c.attribute)
	  & ~MatrixType::Valid) == 0; }

/*************** Make a copy of a matrix on the heap *********************/

GeneralMatrix* Matrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new Matrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* SymmetricMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new SymmetricMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* UpperTriangularMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new UpperTriangularMatrix(*this);
   MatrixErrorNoSpace(gm); return gm;
}

GeneralMatrix* LowerTriangularMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new LowerTriangularMatrix(*this);
   MatrixErrorNoSpace(gm); return gm;
}

GeneralMatrix* DiagonalMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new DiagonalMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* RowVector::Image() const
{
   REPORT
   GeneralMatrix* gm = new RowVector(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* ColumnVector::Image() const
{
   REPORT
   GeneralMatrix* gm = new ColumnVector(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* BandMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new BandMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* UpperBandMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new UpperBandMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* LowerBandMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new LowerBandMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* SymmetricBandMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new SymmetricBandMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* nricMatrix::Image() const
{
   REPORT
   GeneralMatrix* gm = new nricMatrix(*this); MatrixErrorNoSpace(gm);
   return gm;
}

GeneralMatrix* GeneralMatrix::Image() const
{
   REPORT
   Throw(InternalException("Cannot apply Image to this matrix type"));
   return 0;
}


/************************* nricMatrix routines *****************************/

void nricMatrix::MakeRowPointer()
{
   row_pointer = new Real* [nrows]; MatrixErrorNoSpace(row_pointer);
   Real* s = Store() - 1; int i = nrows; Real** rp = row_pointer;
   while (i--) { *rp++ = s; s+=ncols; }
}

void nricMatrix::DeleteRowPointer()
#ifdef Version21
{ if (nrows) delete [] row_pointer; }
#else
{ if (nrows) delete [nrows] row_pointer; }
#endif

void GeneralMatrix::CheckStore() const
{
   if (!store) 
      Throw(ProgramException("NRIC accessing matrix with unset dimensions"));
}


/***************************** CleanUp routines *****************************/

void GeneralMatrix::CleanUp()
{
   // set matrix dimensions to zero, delete storage
   REPORT
   if (store && storage)
   {
      MONITOR_REAL_DELETE("Free (CleanUp)    ",storage,store)
#ifdef Version21
         REPORT delete [] store;
#else
         REPORT delete [storage] store;
#endif
   }
   store=0; storage=0; nrows=0; ncols=0;
}

void nricMatrix::CleanUp()
{ DeleteRowPointer(); GeneralMatrix::CleanUp(); }

void RowVector::CleanUp()
{ GeneralMatrix::CleanUp(); nrows=1; }

void ColumnVector::CleanUp()
{ GeneralMatrix::CleanUp(); ncols=1; }

void CroutMatrix::CleanUp()
{ 
#ifdef Version21
   if (nrows) delete [] indx;
#else
   if (nrows) delete [nrows] indx;
#endif
   GeneralMatrix::CleanUp();
}

void BandLUMatrix::CleanUp()
{ 
#ifdef Version21
   if (nrows) delete [] indx;
   if (storage2) delete [] store2;
#else
   if (nrows) delete [nrows] indx;
   if (storage2) delete [storage2] store2;
#endif
   GeneralMatrix::CleanUp();
}


