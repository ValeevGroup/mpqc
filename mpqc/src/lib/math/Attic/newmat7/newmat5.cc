//$$ newmat5.cxx         Transpose, evaluate etc

// Copyright (C) 1991,2,3: R B Davies

#include "include.h"

#include "newmat.h"
#include "newmatrc.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,5); ExeCount++; }

#define REPORT {}


/************************ carry out operations ******************************/


GeneralMatrix* GeneralMatrix::Transpose(TransposedMatrix* tm, MatrixType mt)
{
   GeneralMatrix* gm1;

   if (Compare(Type().t(),mt))
   {
      REPORT
      gm1 = mt.New(ncols,nrows,tm);
      for (int i=0; i<ncols; i++)
      {
	 MatrixRow mr(gm1, StoreOnExit+DirectPart, i);
         MatrixCol mc(this, mr.Store(), LoadOnEntry, i);
      }
   }
   else
   {
      REPORT
      gm1 = mt.New(ncols,nrows,tm);
      MatrixRow mr(this, LoadOnEntry);
      MatrixCol mc(gm1, StoreOnExit+DirectPart);
      int i = nrows;
      while (i--) { mc.Copy(mr); mr.Next(); mc.Next(); }
   }
   tDelete(); gm1->ReleaseAndDelete(); return gm1;
}

GeneralMatrix* SymmetricMatrix::Transpose(TransposedMatrix*, MatrixType mt)
{ REPORT  return Evaluate(mt); }


GeneralMatrix* DiagonalMatrix::Transpose(TransposedMatrix*, MatrixType mt)
{ REPORT return Evaluate(mt); }

Boolean GeneralMatrix::IsZero() const
{
   REPORT
   Real* s=store; int i=storage;
   while (i--) { if (*s++) return FALSE; }
   return TRUE;
}

GeneralMatrix* ColumnVector::Transpose(TransposedMatrix*, MatrixType mt)
{
   REPORT
   GeneralMatrix* gmx = new RowVector; MatrixErrorNoSpace(gmx);
   gmx->nrows = 1; gmx->ncols = gmx->storage = storage;
   return BorrowStore(gmx,mt);
}

GeneralMatrix* RowVector::Transpose(TransposedMatrix*, MatrixType mt)
{
   REPORT
   GeneralMatrix* gmx = new ColumnVector; MatrixErrorNoSpace(gmx);
   gmx->ncols = 1; gmx->nrows = gmx->storage = storage;
   return BorrowStore(gmx,mt);
}

GeneralMatrix* GeneralMatrix::Evaluate(MatrixType mt)
{
   if (Compare(this->Type(),mt)) { REPORT return this; }
   REPORT
   GeneralMatrix* gmx = mt.New(nrows,ncols,this);
   MatrixRow mr(this, LoadOnEntry);
   MatrixRow mrx(gmx, StoreOnExit+DirectPart);
   int i=nrows;
   while (i--) { mrx.Copy(mr); mrx.Next(); mr.Next(); }
   tDelete(); gmx->ReleaseAndDelete(); return gmx;
}

GeneralMatrix* ConstMatrix::Evaluate(MatrixType mt)
{
   if (Compare(cgm->Type(),mt))
   {
      REPORT
#ifdef TEMPS_DESTROYED_QUICKLY
      GeneralMatrix* gmx = (GeneralMatrix*)cgm; delete this; return gmx;
#else
      return (GeneralMatrix*)cgm;
#endif
   }
   REPORT
   GeneralMatrix* gmx = cgm->Type().New(cgm->Nrows(),cgm->Ncols(),this);
   MatrixRow mr((GeneralMatrix*)cgm, LoadOnEntry);//assume won't change this
   MatrixRow mrx(gmx, StoreOnExit+DirectPart);
   int i=cgm->Nrows();
   while (i--) { mrx.Copy(mr); mrx.Next(); mr.Next(); }
   gmx->ReleaseAndDelete();
#ifdef TEMPS_DESTROYED_QUICKLY
   delete this;
#endif
   return gmx; // no tDelete
}

GeneralMatrix* ShiftedMatrix::Evaluate(MatrixType mt)
{
   gm=((BaseMatrix*&)bm)->Evaluate();
   int nr=gm->Nrows(); int nc=gm->Ncols();
   Compare(gm->Type().AddEqualEl(),mt);
   if (!(mt==gm->Type()))
   {
      REPORT
      GeneralMatrix* gmx = mt.New(nr,nc,this);
      MatrixRow mr(gm, LoadOnEntry); 
      MatrixRow mrx(gmx, StoreOnExit+DirectPart);
      while (nr--) { mrx.Add(mr,f); mrx.Next(); mr.Next(); }
      gmx->ReleaseAndDelete(); gm->tDelete();
#ifdef TEMPS_DESTROYED_QUICKLY
      delete this;
#endif
      return gmx;
   }
   else if (gm->reuse())
   {
      REPORT gm->Add(f);
#ifdef TEMPS_DESTROYED_QUICKLY
      GeneralMatrix* gmx = gm; delete this; return gmx;
#else
      return gm;
#endif
   }
   else
   {
      REPORT GeneralMatrix* gmy = gm->Type().New(nr,nc,this);
      gmy->ReleaseAndDelete(); gmy->Add(gm,f);
#ifdef TEMPS_DESTROYED_QUICKLY
      delete this;
#endif
      return gmy;
   }
}

GeneralMatrix* ScaledMatrix::Evaluate(MatrixType mt)
{
   gm=((BaseMatrix*&)bm)->Evaluate();
   int nr=gm->Nrows(); int nc=gm->Ncols();
   if (Compare(gm->Type(),mt))
   {
      if (gm->reuse())
      {
         REPORT gm->Multiply(f);
#ifdef TEMPS_DESTROYED_QUICKLY
         GeneralMatrix* gmx = gm; delete this; return gmx;
#else
         return gm;
#endif
      }
      else
      {
         REPORT GeneralMatrix* gmx = gm->Type().New(nr,nc,this);
         gmx->ReleaseAndDelete(); gmx->Multiply(gm,f);
#ifdef TEMPS_DESTROYED_QUICKLY
         delete this;
#endif
         return gmx;
      }
   }
   else
   {
      REPORT
      GeneralMatrix* gmx = mt.New(nr,nc,this);
      MatrixRow mr(gm, LoadOnEntry); 
      MatrixRow mrx(gmx, StoreOnExit+DirectPart);
      while (nr--) { mrx.Multiply(mr,f); mrx.Next(); mr.Next(); }
      gmx->ReleaseAndDelete(); gm->tDelete();
#ifdef TEMPS_DESTROYED_QUICKLY
      delete this;
#endif
      return gmx;
   }
}

GeneralMatrix* NegatedMatrix::Evaluate(MatrixType mt)
{
   gm=((BaseMatrix*&)bm)->Evaluate();
   int nr=gm->Nrows(); int nc=gm->Ncols();
   if (Compare(gm->Type(),mt))
   {
      if (gm->reuse())
      {
         REPORT gm->Negate();
#ifdef TEMPS_DESTROYED_QUICKLY
         GeneralMatrix* gmx = gm; delete this; return gmx;
#else
         return gm;
#endif
      }
      else
      {
         REPORT
         GeneralMatrix* gmx = gm->Type().New(nr,nc,this);
         gmx->ReleaseAndDelete(); gmx->Negate(gm);
#ifdef TEMPS_DESTROYED_QUICKLY
         delete this;
#endif
         return gmx;
      }
   }
   else
   {
      REPORT
      GeneralMatrix* gmx = mt.New(nr,nc,this);
      MatrixRow mr(gm, LoadOnEntry); 
      MatrixRow mrx(gmx, StoreOnExit+DirectPart);
      while (nr--) { mrx.Negate(mr); mrx.Next(); mr.Next(); }
      gmx->ReleaseAndDelete(); gm->tDelete();
#ifdef TEMPS_DESTROYED_QUICKLY
      delete this;
#endif
      return gmx;
   }
}   

GeneralMatrix* TransposedMatrix::Evaluate(MatrixType mt)
{
   REPORT
   gm=((BaseMatrix*&)bm)->Evaluate();
   Compare(gm->Type().t(),mt);
   GeneralMatrix* gmx=gm->Transpose(this, mt);
#ifdef TEMPS_DESTROYED_QUICKLY
   delete this;
#endif
   return gmx;
}
   
GeneralMatrix* RowedMatrix::Evaluate(MatrixType mt)
{
   gm = ((BaseMatrix*&)bm)->Evaluate();
   GeneralMatrix* gmx = new RowVector; MatrixErrorNoSpace(gmx);
   gmx->nrows = 1; gmx->ncols = gmx->storage = gm->storage;
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmy = gm; delete this; return gmy->BorrowStore(gmx,mt);
#else
   return gm->BorrowStore(gmx,mt);
#endif
}

GeneralMatrix* ColedMatrix::Evaluate(MatrixType mt)
{
   gm = ((BaseMatrix*&)bm)->Evaluate();
   GeneralMatrix* gmx = new ColumnVector; MatrixErrorNoSpace(gmx);
   gmx->ncols = 1; gmx->nrows = gmx->storage = gm->storage;
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmy = gm; delete this; return gmy->BorrowStore(gmx,mt);
#else
   return gm->BorrowStore(gmx,mt);
#endif
}

GeneralMatrix* DiagedMatrix::Evaluate(MatrixType mt)
{
   gm = ((BaseMatrix*&)bm)->Evaluate();
   GeneralMatrix* gmx = new DiagonalMatrix; MatrixErrorNoSpace(gmx);
   gmx->nrows = gmx->ncols = gmx->storage = gm->storage;
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmy = gm; delete this; return gmy->BorrowStore(gmx,mt);
#else
   return gm->BorrowStore(gmx,mt);
#endif
}

GeneralMatrix* MatedMatrix::Evaluate(MatrixType mt)
{
   Tracer tr("MatedMatrix::Evaluate");
   gm = ((BaseMatrix*&)bm)->Evaluate();
   GeneralMatrix* gmx = new Matrix; MatrixErrorNoSpace(gmx);
   gmx->nrows = nr; gmx->ncols = nc; gmx->storage = gm->storage;
   if (nr*nc != gmx->storage)
      Throw(IncompatibleDimensionsException());
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmy = gm; delete this; return gmy->BorrowStore(gmx,mt);
#else
   return gm->BorrowStore(gmx,mt);
#endif
}

GeneralMatrix* GetSubMatrix::Evaluate(MatrixType mt)
{
   REPORT
   Tracer tr("SubMatrix(evaluate)");
   gm = ((BaseMatrix*&)bm)->Evaluate();
   if (row_number < 0) row_number = gm->Nrows();
   if (col_number < 0) col_number = gm->Ncols();
   if (row_skip+row_number > gm->Nrows() || col_skip+col_number > gm->Ncols())
      Throw(SubMatrixDimensionException());
   if (IsSym) Compare(gm->Type().ssub(), mt); 
   else Compare(gm->Type().sub(), mt);
   GeneralMatrix* gmx = mt.New(row_number, col_number, this);
   int i = row_number;
   MatrixRow mr(gm, LoadOnEntry, row_skip); 
   MatrixRow mrx(gmx, StoreOnExit+DirectPart);
   MatrixRowCol sub;
   while (i--)
   {
      mr.SubRowCol(sub, col_skip, col_number);   // put values in sub
      mrx.Copy(sub); mrx.Next(); mr.Next();
   }
   gmx->ReleaseAndDelete(); gm->tDelete();
#ifdef TEMPS_DESTROYED_QUICKLY
   delete this;
#endif
   return gmx;
}   


GeneralMatrix* ReturnMatrixX::Evaluate(MatrixType mt)
{
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmx = gm; delete this; return gmx->Evaluate(mt);
#else
   return gm->Evaluate(mt);
#endif
}


void GeneralMatrix::Add(GeneralMatrix* gm1, Real f)
{
   REPORT
   Real* s1=gm1->store; Real* s=store; int i=(storage >> 2);
   while (i--)
   { *s++ = *s1++ + f; *s++ = *s1++ + f; *s++ = *s1++ + f; *s++ = *s1++ + f; }
   i = storage & 3; while (i--) *s++ = *s1++ + f;
}
   
void GeneralMatrix::Add(Real f)
{
   REPORT
   Real* s=store; int i=(storage >> 2);
   while (i--) { *s++ += f; *s++ += f; *s++ += f; *s++ += f; }
   i = storage & 3; while (i--) *s++ += f;
}
   
void GeneralMatrix::Negate(GeneralMatrix* gm1)
{
   // change sign of elements
   REPORT
   Real* s1=gm1->store; Real* s=store; int i=(storage >> 2);
   while (i--)
   { *s++ = -(*s1++); *s++ = -(*s1++); *s++ = -(*s1++); *s++ = -(*s1++); }
   i = storage & 3; while(i--) *s++ = -(*s1++);
}
   
void GeneralMatrix::Negate()
{
   REPORT
   Real* s=store; int i=(storage >> 2);
   while (i--)
   { *s = -(*s); s++; *s = -(*s); s++; *s = -(*s); s++; *s = -(*s); s++; }
   i = storage & 3; while(i--) { *s = -(*s); s++; }
}
   
void GeneralMatrix::Multiply(GeneralMatrix* gm1, Real f)
{
   REPORT
   Real* s1=gm1->store; Real* s=store;  int i=(storage >> 2);
   while (i--)
   { *s++ = *s1++ * f; *s++ = *s1++ * f; *s++ = *s1++ * f; *s++ = *s1++ * f; }
   i = storage & 3; while (i--) *s++ = *s1++ * f;
}
   
void GeneralMatrix::Multiply(Real f)
{
   REPORT
   Real* s=store; int i=(storage >> 2);
   while (i--) { *s++ *= f; *s++ *= f; *s++ *= f; *s++ *= f; }
   i = storage & 3; while (i--) *s++ *= f;
}
   

/************************ MatrixInput routines ****************************/

int MatrixInput::n = 0;       // number values still to be read
Real* MatrixInput::r;         // pointer to last thing read
int MatrixInput::depth = 0;   // number of objects of this class in existence

MatrixInput MatrixInput::operator<<(Real f)
{
   if (!(n--))
   { depth=0;  n=0; Throw(ProgramException("List of values too long")); }
   *(++r) = f;
   return MatrixInput();
}


MatrixInput BandMatrix::operator<<(Real)
{
   Throw(ProgramException("Cannot use list read with a BandMatrix"));
   return MatrixInput(); 
}

void BandMatrix::operator<<(const Real*)
{ Throw(ProgramException("Cannot use array read with a BandMatrix")); }

MatrixInput GeneralMatrix::operator<<(Real f)
{
   if (MatrixInput::n)
   {
      MatrixInput::depth=0;            // so we can recover
      MatrixInput::n=0;                // so we can recover
      Throw(ProgramException("A list of values was too short"));
   }
   MatrixInput::n = Storage();
   if (MatrixInput::n<=0)
      Throw(ProgramException("Loading data to zero dimension matrix"));
   MatrixInput::r = Store(); *(MatrixInput::r) = f; MatrixInput::n--;
   return MatrixInput(); 
}

MatrixInput::~MatrixInput()
{
   if (depth && --depth == 0 && n != 0)
   {
      depth = 0; n = 0;                // so we can recover
      Throw(ProgramException("A list of values was too short"));
   }
}

