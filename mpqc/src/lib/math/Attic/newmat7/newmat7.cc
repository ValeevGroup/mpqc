//$$ newmat7.cxx     Invert, solve, binary operations

// Copyright (C) 1991,2,3: R B Davies

#include "include.h"

#include "newmat.h"
#include "newmatrc.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,7); ExeCount++; }

#define REPORT {}


/***************************** solve routines ******************************/

GeneralMatrix* GeneralMatrix::MakeSolver()
{
   REPORT
   GeneralMatrix* gm = new CroutMatrix(*this);
   MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
}

GeneralMatrix* Matrix::MakeSolver()
{
   REPORT
   GeneralMatrix* gm = new CroutMatrix(*this);
   MatrixErrorNoSpace(gm); gm->ReleaseAndDelete(); return gm;
}

void CroutMatrix::Solver(MatrixRowCol& mcout, const MatrixRowCol& mcin)
{
   REPORT
   Real* el = mcin.store; int i = mcin.skip;
   while (i--) *el++ = 0.0;
   el += mcin.storage; i = nrows - mcin.skip - mcin.storage;
   while (i--) *el++ = 0.0;
   lubksb(mcin.store, mcout.skip);
}


// Do we need check for entirely zero output?

void UpperTriangularMatrix::Solver(MatrixRowCol& mcout,
   const MatrixRowCol& mcin)
{
   REPORT
   Real* elx = mcin.store+mcout.skip; int i = mcin.skip-mcout.skip;
   while (i-- > 0) *elx++ = 0.0;
   int nr = mcin.skip+mcin.storage; elx = mcin.store+nr; Real* el = elx;
   int j = mcout.skip+mcout.storage-nr; int nc = ncols-nr; i = nr-mcout.skip;
   while (j-- > 0) *elx++ = 0.0;
   Real* Ael = store + (nr*(2*ncols-nr+1))/2; j = 0;
   while (i-- > 0)
   {
      elx = el; Real sum = 0.0; int jx = j++; Ael -= nc;
      while (jx--) sum += *(--Ael) * *(--elx);
      elx--; *elx = (*elx - sum) / *(--Ael);
   }
}

void LowerTriangularMatrix::Solver(MatrixRowCol& mcout,
   const MatrixRowCol& mcin)
{
   REPORT
   Real* elx = mcin.store+mcout.skip; int i = mcin.skip-mcout.skip;
   while (i-- > 0) *elx++ = 0.0;
   int nc = mcin.skip; i = nc+mcin.storage; elx = mcin.store+i;
   int nr = mcout.skip+mcout.storage; int j = nr-i; i = nr-nc;
   while (j-- > 0) *elx++ = 0.0;
   Real* el = mcin.store+nc; Real* Ael = store + (nc*(nc+1))/2; j = 0;
   while (i-- > 0)
   {
      elx = el; Real sum = 0.0; int jx = j++; Ael += nc;
      while (jx--) sum += *Ael++ * *elx++;
      *elx = (*elx - sum) / *Ael++;
   }
}

/******************* carry out binary operations *************************/

static GeneralMatrix*
   GeneralAdd(GeneralMatrix*,GeneralMatrix*,AddedMatrix*,MatrixType);
static GeneralMatrix*
   GeneralSub(GeneralMatrix*,GeneralMatrix*,SubtractedMatrix*,MatrixType);
static GeneralMatrix*
   GeneralMult(GeneralMatrix*,GeneralMatrix*,MultipliedMatrix*,MatrixType);
static GeneralMatrix*
   GeneralSolv(GeneralMatrix*,GeneralMatrix*,BaseMatrix*,MatrixType);

GeneralMatrix* AddedMatrix::Evaluate(MatrixType mt)
{
   REPORT
   gm1=((BaseMatrix*&)bm1)->Evaluate();
   gm2=((BaseMatrix*&)bm2)->Evaluate();
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmx;
   Try { gmx = GeneralAdd(gm1,gm2,this,mt); }
   CatchAll { delete this; Throw(); }
   delete this; return gmx;
#else
   return GeneralAdd(gm1,gm2,this,mt);
#endif   
}

GeneralMatrix* SubtractedMatrix::Evaluate(MatrixType mt)
{
   REPORT
   gm1=((BaseMatrix*&)bm1)->Evaluate();
   gm2=((BaseMatrix*&)bm2)->Evaluate();
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmx;
   Try { gmx = GeneralSub(gm1,gm2,this,mt); }
   CatchAll { delete this; Throw(); }
   delete this; return gmx;
#else
   return GeneralSub(gm1,gm2,this,mt);
#endif   
}

GeneralMatrix* MultipliedMatrix::Evaluate(MatrixType mt)
{
   REPORT
   gm2 = ((BaseMatrix*&)bm2)->Evaluate();
   gm2 = gm2->Evaluate(gm2->Type().MultRHS());     // no symmetric on RHS
   gm1=((BaseMatrix*&)bm1)->Evaluate();
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmx;
   Try { gmx = GeneralMult(gm1, gm2, this, mt); }
   CatchAll { delete this; Throw(); }
   delete this; return gmx;
#else
   return GeneralMult(gm1, gm2, this, mt);
#endif   
}

GeneralMatrix* SolvedMatrix::Evaluate(MatrixType mt)
{
   REPORT
   gm1=((BaseMatrix*&)bm1)->Evaluate();
   gm2=((BaseMatrix*&)bm2)->Evaluate();
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmx;
   Try { gmx = GeneralSolv(gm1,gm2,this,mt); }
   CatchAll { delete this; Throw(); }
   delete this; return gmx;
#else
   return GeneralSolv(gm1,gm2,this,mt);
#endif   
}

// routines for adding or subtracting matrices of identical storage structure

static void Add(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
{
   REPORT
   Real* s1=gm1->Store(); Real* s2=gm2->Store();
   Real* s=gm->Store(); int i=gm->Storage() >> 2;
   while (i--)
   {
       *s++ = *s1++ + *s2++; *s++ = *s1++ + *s2++;
       *s++ = *s1++ + *s2++; *s++ = *s1++ + *s2++;
   }
   i=gm->Storage() & 3; while (i--) *s++ = *s1++ + *s2++;
}
   
static void Add(GeneralMatrix* gm, GeneralMatrix* gm2)
{
   REPORT
   Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
   while (i--)
   { *s++ += *s2++; *s++ += *s2++; *s++ += *s2++; *s++ += *s2++; }
   i=gm->Storage() & 3; while (i--) *s++ += *s2++;
}

static void Subtract(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
{
   REPORT
   Real* s1=gm1->Store(); Real* s2=gm2->Store();
   Real* s=gm->Store(); int i=gm->Storage() >> 2;
   while (i--)
   {
       *s++ = *s1++ - *s2++; *s++ = *s1++ - *s2++;
       *s++ = *s1++ - *s2++; *s++ = *s1++ - *s2++;
   }
   i=gm->Storage() & 3; while (i--) *s++ = *s1++ - *s2++;
}

static void Subtract(GeneralMatrix* gm, GeneralMatrix* gm2)
{
   REPORT
   Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
   while (i--)
   { *s++ -= *s2++; *s++ -= *s2++; *s++ -= *s2++; *s++ -= *s2++; }
   i=gm->Storage() & 3; while (i--) *s++ -= *s2++;
}

static void ReverseSubtract(GeneralMatrix* gm, GeneralMatrix* gm2)
{
   REPORT
   Real* s2=gm2->Store(); Real* s=gm->Store(); int i=gm->Storage() >> 2;
   while (i--)
   {
      *s = *s2++ - *s; s++; *s = *s2++ - *s; s++;
      *s = *s2++ - *s; s++; *s = *s2++ - *s; s++;
   }
   i=gm->Storage() & 3; while (i--) { *s = *s2++ - *s; s++; }
}

// routines for adding or subtracting matrices of different storage structure

static void AddDS(GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
{
   int nr = gm->Nrows();
   MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
   MatrixRow mr(gm, StoreOnExit+DirectPart);
   while (nr--) { mr.Add(mr1,mr2); mr1.Next(); mr2.Next(); mr.Next(); }
}
   
static void AddDS(GeneralMatrix* gm, GeneralMatrix* gm2)
// Add into first argument
{
   int nr = gm->Nrows();
   MatrixRow mr(gm, StoreOnExit+LoadOnEntry+DirectPart);
   MatrixRow mr2(gm2, LoadOnEntry);
   while (nr--) { mr.Add(mr2); mr.Next(); mr2.Next(); }
}

static void SubtractDS
   (GeneralMatrix* gm, GeneralMatrix* gm1, GeneralMatrix* gm2)
{
   int nr = gm->Nrows();
   MatrixRow mr1(gm1, LoadOnEntry); MatrixRow mr2(gm2, LoadOnEntry);
   MatrixRow mr(gm, StoreOnExit+DirectPart);
   while (nr--) { mr.Sub(mr1,mr2); mr1.Next(); mr2.Next(); mr.Next(); }
}

static void SubtractDS(GeneralMatrix* gm, GeneralMatrix* gm2)
{
   int nr = gm->Nrows();
   MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart);
   MatrixRow mr2(gm2, LoadOnEntry);
   while (nr--) { mr.Sub(mr2); mr.Next(); mr2.Next(); }
}

static void ReverseSubtractDS(GeneralMatrix* gm, GeneralMatrix* gm2)
{
   int nr = gm->Nrows();
   MatrixRow mr(gm, LoadOnEntry+StoreOnExit+DirectPart);
   MatrixRow mr2(gm2, LoadOnEntry);
   while (nr--) { mr.RevSub(mr2); mr2.Next(); mr.Next(); }
}

#ifdef __GNUG__
void AddedMatrix::SelectVersion
   (MatrixType mtx, int& c1, int& c2) const
#else
void AddedMatrix::SelectVersion
   (MatrixType mtx, Boolean& c1, Boolean& c2) const
#endif
// for determining version of add and subtract routines
// will need to modify if further matrix structures are introduced
{
   MatrixBandWidth bm1 = gm1->BandWidth();
   MatrixBandWidth bm2 = gm2->BandWidth();
   MatrixBandWidth bmx = bm1 + bm2;
   c1 = (mtx == gm1->Type()) && (bmx == bm1);
   c2 = (mtx == gm2->Type()) && (bmx == bm2);
}

static GeneralMatrix* GeneralAdd(GeneralMatrix* gm1, GeneralMatrix* gm2,
   AddedMatrix* am, MatrixType mtx)
{
   Tracer tr("GeneralAdd");
   int nr=gm1->Nrows(); int nc=gm1->Ncols();
   if (nr!=gm2->Nrows() || nc!=gm2->Ncols()) 
      Throw(IncompatibleDimensionsException());
   Compare(gm1->Type() + gm2->Type(),mtx);
#ifdef __GNUG__
   int c1,c2; am->SelectVersion(mtx,c1,c2);
#else
   Boolean c1,c2; am->SelectVersion(mtx,c1,c2); // causes problems for g++
#endif
   if (c1 && c2)
   {
      if (gm1->reuse()) { REPORT Add(gm1,gm2); gm2->tDelete(); return gm1; }
      else if (gm2->reuse()) { REPORT Add(gm2,gm1); return gm2; }
      else
      {
         REPORT GeneralMatrix* gmx = gm1->Type().New(nr,nc,am);
         gmx->ReleaseAndDelete(); Add(gmx,gm1,gm2); return gmx;
      }
   }
   else
   {
      if (c1 && gm1->reuse() )               // must have type test first
      { REPORT AddDS(gm1,gm2); gm2->tDelete(); return gm1; }
      else if (c2 && gm2->reuse() )
      { REPORT AddDS(gm2,gm1); if (!c1) gm1->tDelete(); return gm2; }
      else
      {
         REPORT
	 GeneralMatrix* gmx = mtx.New(nr,nc,am); AddDS(gmx,gm1,gm2);
	 if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
         gmx->ReleaseAndDelete(); return gmx;
      }
   }
}


static GeneralMatrix* GeneralSub(GeneralMatrix* gm1, GeneralMatrix* gm2,
   SubtractedMatrix* sm, MatrixType mtx)
{
   Tracer tr("GeneralSub");
   Compare(gm1->Type() + gm2->Type(),mtx);
   int nr=gm1->Nrows(); int nc=gm1->Ncols();
   if (nr!=gm2->Nrows() || nc!=gm2->Ncols())
      Throw(IncompatibleDimensionsException());
#ifdef __GNUG__
   int c1,c2; sm->SelectVersion(mtx,c1,c2);
#else
   Boolean c1,c2; sm->SelectVersion(mtx,c1,c2); // causes problems for g++
#endif
   if (c1 && c2)
   {
      if (gm1->reuse())
      { REPORT Subtract(gm1,gm2); gm2->tDelete(); return gm1; }
      else if (gm2->reuse()) { REPORT ReverseSubtract(gm2,gm1); return gm2; }
      else
      {
         REPORT
	 GeneralMatrix* gmx = gm1->Type().New(nr,nc,sm);
         gmx->ReleaseAndDelete(); Subtract(gmx,gm1,gm2); return gmx;
      }
   }
   else
   {
      if ( c1 && gm1->reuse() )
      { REPORT  SubtractDS(gm1,gm2); gm2->tDelete(); return gm1; }
      else if ( c2 && gm2->reuse() )
      {
         REPORT
         ReverseSubtractDS(gm2,gm1); if (!c1) gm1->tDelete(); return gm2;
      }
      else
      {
         REPORT
	 GeneralMatrix* gmx = mtx.New(nr,nc,sm); SubtractDS(gmx,gm1,gm2);
	 if (!c1) gm1->tDelete(); if (!c2) gm2->tDelete();
	 gmx->ReleaseAndDelete(); return gmx;
      }
   }
}

static GeneralMatrix* GeneralMult1(GeneralMatrix* gm1, GeneralMatrix* gm2,
   MultipliedMatrix* mm, MatrixType mtx)
{
   REPORT
   Tracer tr("GeneralMult1");
   int nr=gm1->Nrows(); int nc=gm2->Ncols();
   if (gm1->Ncols() !=gm2->Nrows())
      Throw(IncompatibleDimensionsException());
   GeneralMatrix* gmx = mtx.New(nr,nc,mm);

   MatrixCol mcx(gmx, StoreOnExit+DirectPart);
   MatrixCol mc2(gm2, LoadOnEntry);
   while (nc--)
   {
      MatrixRow mr1(gm1, LoadOnEntry, mcx.Skip());
      Real* el = mcx();                              // pointer to an element
      int n = mcx.Storage();
      while (n--) { *(el++) = DotProd(mr1,mc2); mr1.Next(); }
      mc2.Next(); mcx.Next();
   }
   gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
}

static GeneralMatrix* GeneralMult2(GeneralMatrix* gm1, GeneralMatrix* gm2,
   MultipliedMatrix* mm, MatrixType mtx)
{
   // version that accesses by row only - not good for thin matrices
   // or column vectors in right hand term. Needs fixing
   REPORT
   Tracer tr("GeneralMult2");
   int nr=gm1->Nrows(); int nc=gm2->Ncols();
   if (gm1->Ncols() !=gm2->Nrows())
      Throw(IncompatibleDimensionsException());
   GeneralMatrix* gmx = mtx.New(nr,nc,mm);

   Real* el = gmx->Store(); int n = gmx->Storage();
   while (n--) *el++ = 0.0;
   MatrixRow mrx(gmx, LoadOnEntry+StoreOnExit+DirectPart);
   MatrixRow mr1(gm1, LoadOnEntry);
   while (nr--)
   {
      MatrixRow mr2(gm2, LoadOnEntry, mr1.Skip());
      el = mr1();                              // pointer to an element
      n = mr1.Storage();
      while (n--) { mrx.AddScaled(mr2, *el++); mr2.Next(); }
      mr1.Next(); mrx.Next();
   }
   gmx->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gmx;
}

static GeneralMatrix* mmMult(GeneralMatrix* gm1, GeneralMatrix* gm2)
{
   // matrix multiplication for type Matrix only
   REPORT
   Tracer tr("MatrixMult");

   int nr=gm1->Nrows(); int ncr=gm1->Ncols(); int nc=gm2->Ncols();
   if (ncr != gm2->Nrows()) Throw(IncompatibleDimensionsException());

   Matrix* gm = new Matrix(nr,nc); MatrixErrorNoSpace(gm);

   Real* s1=gm1->Store(); Real* s2=gm2->Store(); Real* s=gm->Store();
   
   if (ncr)
   {
      while (nr--)
      {
         Real* s2x = s2; int j = ncr;
         Real* sx = s; Real f = *s1++; int k = nc;
         while (k--) *sx++ = f * *s2x++;
         while (--j)
            { sx = s; f = *s1++; k = nc; while (k--) *sx++ += f * *s2x++; }
         s = sx;
      }
   }
   else *gm = 0.0;

   gm->ReleaseAndDelete(); gm1->tDelete(); gm2->tDelete(); return gm;
}

static GeneralMatrix* GeneralMult(GeneralMatrix* gm1, GeneralMatrix* gm2,
   MultipliedMatrix* mm, MatrixType mtx)
{
   if ( Rectangular(gm1->Type(), gm2->Type(), mtx)) return mmMult(gm1, gm2);
   else
   {
      Compare(gm1->Type() * gm2->Type(),mtx);
      int nr = gm2->Nrows(); int nc = gm2->Ncols();
      if (nc <= 5 && nr > nc) return GeneralMult1(gm1, gm2, mm, mtx);
      else return GeneralMult2(gm1, gm2, mm, mtx);
   }
}

static GeneralMatrix* GeneralSolv(GeneralMatrix* gm1, GeneralMatrix* gm2,
   BaseMatrix* sm, MatrixType mtx)
{
   REPORT
   Tracer tr("GeneralSolv");
   Compare(gm1->Type().i() * gm2->Type(),mtx);
   int nr=gm1->Nrows(); int nc=gm2->Ncols();
   if (gm1->Ncols() !=gm2->Nrows())
      Throw(IncompatibleDimensionsException());
   GeneralMatrix* gmx = mtx.New(nr,nc,sm);

   Real* r = new Real [nr]; MatrixErrorNoSpace(r);
#ifndef ATandT
   MONITOR_REAL_NEW("Make   (GenSolv)",nr,r)
                           // deleted for ATandT, to balance deletion below
#endif
   MatrixCol mcx(gmx, r, StoreOnExit+DirectPart);   // copy to and from r
   MatrixCol mc2(gm2, r, LoadOnEntry);
   GeneralMatrix* gms = gm1->MakeSolver();
   Try
   {
      while (nc--) { gms->Solver(mcx, mc2); mcx.Next(); mc2.Next(); }
   }
   CatchAll
   {
      gms->tDelete(); delete gmx; gm2->tDelete();
#ifndef ATandT
      MONITOR_REAL_DELETE("Delete (GenSolv)",nr,r)
                          // ATandT version 2.1 gives an internal error
#endif
#ifdef Version21
      delete [] r;
#else
      delete [nr] r;
#endif
      Throw();
   }
   gms->tDelete(); gmx->ReleaseAndDelete(); gm2->tDelete();
#ifndef ATandT
   MONITOR_REAL_DELETE("Delete (GenSolv)",nr,r)
                          // ATandT version 2.1 gives an internal error
#endif
#ifdef Version21
   delete [] r;
#else
   delete [nr] r;
#endif
   return gmx;
}

GeneralMatrix* InvertedMatrix::Evaluate(MatrixType mtx)
{
   // Matrix Inversion - use solve routines
   REPORT
   gm=((BaseMatrix*&)bm)->Evaluate();
   int n = gm->Nrows(); DiagonalMatrix I(n); I=1.0;
#ifdef TEMPS_DESTROYED_QUICKLY
   GeneralMatrix* gmx;
   Try { gmx = GeneralSolv(gm,&I,this,mtx); }
   CatchAll { delete this; Throw(); }
   delete this; return gmx;
#else
   return GeneralSolv(gm,&I,this,mtx);
#endif   
}


/*************************** norm functions ****************************/

Real BaseMatrix::Norm1() const
{
   // maximum of sum of absolute values of a column
   REPORT
   GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   int nc = gm->Ncols(); Real value = 0.0;
   MatrixCol mc(gm, LoadOnEntry);
   while (nc--)
      { Real v = mc.SumAbsoluteValue(); if (value < v) value = v; mc.Next(); }
   gm->tDelete(); return value;
}

Real BaseMatrix::NormInfinity() const
{
   // maximum of sum of absolute values of a row
   REPORT
   GeneralMatrix* gm = ((BaseMatrix&)*this).Evaluate();
   int nr = gm->Nrows(); Real value = 0.0;
   MatrixRow mr(gm, LoadOnEntry);
   while (nr--)
      { Real v = mr.SumAbsoluteValue(); if (value < v) value = v; mr.Next(); }
   gm->tDelete(); return value;
}
