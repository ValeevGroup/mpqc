//$$ newmat3.cxx        Matrix get and restore rows and columns

// Copyright (C) 1991,2,3: R B Davies


#include "include.h"

#include "newmat.h"
#include "newmatrc.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,3); ExeCount++; }

#define REPORT {}

//#define MONITOR(what,storage,store) \
//   { cout << what << " " << storage << " at " << (long)store << "\n"; }

#define MONITOR(what,store,storage) {}


// Control bits codes for GetRow, GetCol, RestoreRow, RestoreCol
//
// LoadOnEntry:
//    Load data into MatrixRow or Col dummy array under GetRow or GetCol
// StoreOnExit:
//    Restore data to original matrix under RestoreRow or RestoreCol
// IsACopy:
//    Set by GetRow/Col: MatrixRow or Col array is a copy
// DirectPart:
//    Load or restore only part directly stored; must be set with StoreOnExit
//    Still have decide  how to handle this with symmetric
// StoreHere:
//    used in columns only - store data at supplied storage address, adjusted
//    for skip; used for GetCol, NextCol & RestoreCol. No need to fill out
//    zeros.


// These will work as a default
// but need to override NextRow for efficiency

// Assume pointer arithmetic works for pointers out of range - not strict C++.


void GeneralMatrix::NextRow(MatrixRowCol& mrc)
{
   REPORT
   if (+(mrc.cw*StoreOnExit)) { REPORT this->RestoreRow(mrc); }
   if (+(mrc.cw*IsACopy))
   {
      REPORT
      Real* s = mrc.store + mrc.skip;
      MONITOR_REAL_DELETE("Free   (NextRow)",mrc.storage,s)
#ifdef Version21
      delete [] s;
#else
      delete [mrc.storage] s;
#endif
   }
   mrc.rowcol++;
   if (mrc.rowcol<nrows) { REPORT this->GetRow(mrc); }
   else { REPORT mrc.cw -= (StoreOnExit+IsACopy); }
}

void GeneralMatrix::NextCol(MatrixRowCol& mrc)
{
   REPORT                                                // 423
   if (+(mrc.cw*StoreOnExit)) { REPORT this->RestoreCol(mrc); }
   int t1 = +(mrc.cw*IsACopy); int t2 = !(mrc.cw*StoreHere);
   if ( t1 && t2 )
   {
      REPORT                                             // not accessed
      Real* s = mrc.store + mrc.skip;
      MONITOR_REAL_DELETE("Free   (NextCol)",mrc.storage,s) 
#ifdef Version21
      delete [] s;
#else
      delete [mrc.storage] s;
#endif
   }
   mrc.rowcol++;
   if (mrc.rowcol<ncols) { REPORT this->GetCol(mrc); }
   else { REPORT mrc.cw -= (StoreOnExit+IsACopy); }
}


// routines for matrix

void Matrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; mrc.cw-=IsACopy; mrc.storage=ncols;
   mrc.store=store+mrc.rowcol*ncols;
}


void Matrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; mrc.storage=nrows; int t1 = !(mrc.cw*StoreHere);
   if ( ncols==1 && t1 )
      { REPORT mrc.cw-=IsACopy; mrc.store=store; }           // not accessed
   else
   {
      mrc.cw+=IsACopy; Real* ColCopy;
      if ( !(mrc.cw*StoreHere) )
      {
         REPORT
         ColCopy = new Real [nrows]; MatrixErrorNoSpace(ColCopy);
         MONITOR_REAL_NEW("Make (MatGetCol)",nrows,ColCopy)
         mrc.store = ColCopy;
      }
      else { REPORT ColCopy = mrc.store; }
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+mrc.rowcol; int i=nrows;
         while (i--) { *ColCopy++ = *Mstore; Mstore+=ncols; }
      }
   }
}

void Matrix::RestoreCol(MatrixRowCol& mrc)
{
//  if (mrc.cw*StoreOnExit)
   REPORT                                   // 429
   if (+(mrc.cw*IsACopy))
   {
      REPORT                                // 426
      Real* Mstore = store+mrc.rowcol; int i=nrows; Real* Cstore = mrc.store;
      while (i--) { *Mstore = *Cstore++; Mstore+=ncols; }
   }
}

void Matrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrMat(); }  // 1808

void Matrix::NextCol(MatrixRowCol& mrc)
{
   REPORT                                        // 632
   if (+(mrc.cw*StoreOnExit)) { REPORT RestoreCol(mrc); }
   mrc.rowcol++;
   if (mrc.rowcol<ncols)
   {
      if (+(mrc.cw*LoadOnEntry))
      {
	 REPORT
         Real* ColCopy = mrc.store;
         Real* Mstore = store+mrc.rowcol; int i=nrows;
         while (i--) { *ColCopy++ = *Mstore; Mstore+=ncols; }
      }
   }
   else { REPORT mrc.cw -= StoreOnExit; }
}

// routines for diagonal matrix

void DiagonalMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=mrc.rowcol; mrc.cw-=IsACopy; mrc.storage=1; mrc.store=store;
}

void DiagonalMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT 
   mrc.skip=mrc.rowcol; mrc.storage=1;
   if (+(mrc.cw*StoreHere))
      { REPORT *(mrc.store+mrc.rowcol)=*(store+mrc.rowcol); mrc.cw+=IsACopy; }
   else { REPORT mrc.store = store; mrc.cw-=IsACopy; }     // not accessed
}

void DiagonalMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrDiag(); }
						      // 800
void DiagonalMatrix::NextCol(MatrixRowCol& mrc)
{
   REPORT
   if (+(mrc.cw*StoreHere))
   {
      if (+(mrc.cw*StoreOnExit))
         { REPORT *(store+mrc.rowcol)=*(mrc.store+mrc.rowcol); }
      mrc.IncrDiag();
      int t1 = +(mrc.cw*LoadOnEntry);
      if (t1 && mrc.rowcol < ncols)
         { REPORT *(mrc.store+mrc.rowcol)=*(store+mrc.rowcol); }
   }
   else { REPORT mrc.IncrDiag(); }                     // not accessed
}

// routines for upper triangular matrix

void UpperTriangularMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int row = mrc.rowcol; mrc.skip=row; mrc.cw-=IsACopy;
   mrc.storage=ncols-row; mrc.store=store+(row*(2*ncols-row-1))/2;
}


void UpperTriangularMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; mrc.cw+=IsACopy; int i=mrc.rowcol+1; mrc.storage=i;
   Real* ColCopy;
   if ( !(mrc.cw*StoreHere) )
   {
      REPORT                                              // not accessed
      ColCopy = new Real [i]; MatrixErrorNoSpace(ColCopy);
      MONITOR_REAL_NEW("Make (UT GetCol)",i,ColCopy) 
      mrc.store = ColCopy;
   }
   else { REPORT ColCopy = mrc.store; }
   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* Mstore = store+mrc.rowcol; int j = ncols;
      while (i--) { *ColCopy++ = *Mstore; Mstore += --j; }
   }
}

void UpperTriangularMatrix::RestoreCol(MatrixRowCol& mrc)
{
//  if (mrc.cw*StoreOnExit)
  {
     REPORT
     Real* Mstore = store+mrc.rowcol; int i=mrc.rowcol+1; int j = ncols;
     Real* Cstore = mrc.store;
     while (i--) { *Mstore = *Cstore++; Mstore += --j; }
  }
}

void UpperTriangularMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrUT(); }
						      // 722

// routines for lower triangular matrix

void LowerTriangularMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int row=mrc.rowcol; mrc.skip=0; mrc.cw-=IsACopy; mrc.storage=row+1;
   mrc.store=store+(row*(row+1))/2;
}

void LowerTriangularMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   int col=mrc.rowcol; mrc.skip=col; mrc.cw+=IsACopy;
   int i=nrows-col; mrc.storage=i; Real* ColCopy;
   if ( !(mrc.cw*StoreHere) )
   {
      REPORT                                            // not accessed
      ColCopy = new Real [i]; MatrixErrorNoSpace(ColCopy);
      MONITOR_REAL_NEW("Make (LT GetCol)",i,ColCopy) 
      mrc.store = ColCopy-col;
   }
   else { REPORT ColCopy = mrc.store+col; }
   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* Mstore = store+(col*(col+3))/2;
      while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
   }
}

void LowerTriangularMatrix::RestoreCol(MatrixRowCol& mrc)
{
//  if (mrc.cw*StoreOnExit)
   {
      REPORT
      int col=mrc.rowcol; Real* Cstore = mrc.store+col;
      Real* Mstore = store+(col*(col+3))/2; int i=nrows-col;
      while (i--) { *Mstore = *Cstore++; Mstore += ++col; }
   }
}

void LowerTriangularMatrix::NextRow(MatrixRowCol& mrc) { REPORT mrc.IncrLT(); }
					                 //712
// routines for symmetric matrix

void SymmetricMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT                                                //571
   mrc.skip=0; int row=mrc.rowcol;
   if (+(mrc.cw*DirectPart))
   {
      REPORT
      mrc.cw-=IsACopy; mrc.storage=row+1; mrc.store=store+(row*(row+1))/2;
   }
   else
   {
      mrc.cw+=IsACopy; mrc.storage=ncols;
      Real* RowCopy = new Real [ncols]; MatrixErrorNoSpace(RowCopy);
      MONITOR_REAL_NEW("Make (SymGetRow)",ncols,RowCopy) 
      mrc.store = RowCopy;
      if (+(mrc.cw*LoadOnEntry))
      {
	 REPORT                                         // 544
         Real* Mstore = store+(row*(row+1))/2; int i = row;
         while (i--) *RowCopy++ = *Mstore++;
         i = ncols-row;
	 while (i--) { *RowCopy++ = *Mstore; Mstore += ++row; }
      }
   }
}

// need to check this out under StoreHere

void SymmetricMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.skip=0; int col=mrc.rowcol;
   if (+(mrc.cw*DirectPart))
   {
      REPORT                                         // not accessed
      mrc.cw-=IsACopy; mrc.storage=col+1; mrc.store=store+(col*(col+1))/2;
   }
   else
   {
      mrc.cw+=IsACopy; mrc.storage=ncols; Real* ColCopy;
      if ( !(mrc.cw*StoreHere) )
      {
         REPORT                                      // not accessed
         ColCopy = new Real [ncols]; MatrixErrorNoSpace(ColCopy);
         MONITOR_REAL_NEW("Make (SymGetCol)",ncols,ColCopy) 
         mrc.store = ColCopy;
      }
      else { REPORT ColCopy = mrc.store; }
      if (+(mrc.cw*LoadOnEntry))
      {
         REPORT
         Real* Mstore = store+(col*(col+1))/2; int i = col;
         while (i--) *ColCopy++ = *Mstore++;
         i = ncols-col;
	 while (i--) { *ColCopy++ = *Mstore; Mstore += ++col; }
      }
   }
}

//void SymmetricMatrix::RestoreRow(int row, Real* Rstore)
//{
////   if (cw*IsACopy && cw*StoreOnExit)
//   {
//      Real* Mstore = store+(row*(row+1))/2; int i = row+1;
//      while (i--) *Mstore++ = *Rstore++;
//   }
//}

//void SymmetricMatrix::RestoreCol(int col, Real* Cstore)
//{
////   if (cw*IsACopy && cw*StoreOnExit)
//   {
//      Real* Mstore = store+(col*(col+3))/2;
//      int i = nrows-col; int j = col;
//      while (i--) { *Mstore = *Cstore++; Mstore+= ++j; }
//   }
//}

// routines for row vector

void RowVector::GetCol(MatrixRowCol& mrc)
{
   REPORT 
   mrc.skip=0; mrc.storage=1;
   if (mrc.cw >= StoreHere)
   {
      if (mrc.cw >= LoadOnEntry) { REPORT *(mrc.store) = *(store+mrc.rowcol); }
      mrc.cw+=IsACopy;
   }
   else  { REPORT mrc.store = store+mrc.rowcol; mrc.cw-=IsACopy; }
                                                         // not accessed
}

void RowVector::NextCol(MatrixRowCol& mrc) 
{
   REPORT
   if (+(mrc.cw*StoreHere))
   {
      if (+(mrc.cw*StoreOnExit)) { REPORT *(store+mrc.rowcol)=*(mrc.store); }
							 // not accessed
      mrc.rowcol++;
      if (mrc.rowcol < ncols)
      {
	 if (+(mrc.cw*LoadOnEntry)) { REPORT *(mrc.store)=*(store+mrc.rowcol); }
      }
      else { REPORT mrc.cw -= StoreOnExit; }
   }
   else  { REPORT mrc.rowcol++; mrc.store++; }             // not accessed
}

void RowVector::RestoreCol(MatrixRowCol& mrc)
{
   REPORT                                            // not accessed
   if (mrc.cw>=IsACopy)  { REPORT *(store+mrc.rowcol)=*(mrc.store); }
}


// routines for band matrices

void BandMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   mrc.cw -= IsACopy; int r = mrc.rowcol; int w = lower+1+upper;
   int s = r-lower; mrc.store = store+(r*w-s); if (s<0) { w += s; s = 0; }
   mrc.skip = s; s += w-ncols; if (s>0) w -= s; mrc.storage = w;
}

// make special versions of this for upper and lower band matrices

void BandMatrix::NextRow(MatrixRowCol& mrc)
{
   REPORT
   int r = ++mrc.rowcol; mrc.store += lower+upper;
   if (r<=lower) mrc.storage++; else mrc.skip++;
   if (r>=ncols-upper) mrc.storage--;
}

void BandMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   mrc.cw += IsACopy; int c = mrc.rowcol; int n = lower+upper; int w = n+1;
   int b; int s = c-upper; Real* ColCopy;
   if (s<=0) { w += s; s = 0; b = c+lower; } else b = s*w+n;
   mrc.skip = s; s += w-nrows; if (s>0) w -= s; mrc.storage = w;
   if ( !(mrc.cw*StoreHere) )
   {
      REPORT
      ColCopy = new Real [w]; MatrixErrorNoSpace(ColCopy);
      MONITOR_REAL_NEW("Make (BMGetCol)",w,ColCopy)
      mrc.store = ColCopy-mrc.skip;
   }
   else { REPORT ColCopy = mrc.store+mrc.skip; }
   if (+(mrc.cw*LoadOnEntry))
   {
      REPORT
      Real* Mstore = store+b;
      while (w--) { *ColCopy++ = *Mstore; Mstore+=n; }
   }
}

void BandMatrix::RestoreCol(MatrixRowCol& mrc)
{
//  if (mrc.cw*StoreOnExit)
   REPORT
   int c = mrc.rowcol; int n = lower+upper; int s = c-upper;
   Real* Mstore = store + ((s<=0) ? c+lower : s*n+s+n);
   Real* Cstore = mrc.store+mrc.skip; int w = mrc.storage;
   while (w--) { *Mstore = *Cstore++; Mstore += n; }
}

// routines for symmetric band matrix

void SymmetricBandMatrix::GetRow(MatrixRowCol& mrc)
{
   REPORT
   int r=mrc.rowcol; int s = r-lower; int w1 = lower+1; int o = r*w1;
   if (s<0) { w1 += s; o -= s; s = 0; }  mrc.skip = s;

   if (+(mrc.cw*DirectPart))
      { REPORT  mrc.cw -= IsACopy; mrc.store = store+o-s; mrc.storage = w1; }
   else
   {
      mrc.cw += IsACopy; int w = w1+lower; s += w-ncols;
      if (s>0) w -= s; mrc.storage = w; int w2 = w-w1;
      Real* RowCopy = new Real [w]; MatrixErrorNoSpace(RowCopy);
      MONITOR_REAL_NEW("Make (SmBGetRow)",w,RowCopy) 
      mrc.store = RowCopy-mrc.skip;
      if (+(mrc.cw*LoadOnEntry))
      {
	 REPORT
         Real* Mstore = store+o;
         while (w1--) *RowCopy++ = *Mstore++;   Mstore--;
         while (w2--) { Mstore += lower; *RowCopy++ = *Mstore; }
      }
   }
}

// need to check this out under StoreHere

void SymmetricBandMatrix::GetCol(MatrixRowCol& mrc)
{
   REPORT
   int c=mrc.rowcol; int s = c-lower; int w1 = lower+1; int o = c*w1;
   if (s<0) { w1 += s; o -= s; s = 0; }  mrc.skip = s;

   if (+(mrc.cw*DirectPart))
      { REPORT  mrc.cw -= IsACopy; mrc.store = store+o-s; mrc.storage = w1; }
   else
   {
      mrc.cw += IsACopy; int w = w1+lower; s += w-ncols;
      if (s>0) w -= s; mrc.storage = w; int w2 = w-w1; Real* ColCopy;
      if ( !(mrc.cw*StoreHere) )
      {
         ColCopy = new Real [w]; MatrixErrorNoSpace(ColCopy);
         MONITOR_REAL_NEW("Make (SmBGetCol)",w,ColCopy) 
         mrc.store = ColCopy-mrc.skip;
      }
      else { REPORT ColCopy = mrc.store+mrc.skip; }
      if (+(mrc.cw*LoadOnEntry))
      {
	 REPORT
         Real* Mstore = store+o;
         while (w1--) *ColCopy++ = *Mstore++;   Mstore--;
         while (w2--) { Mstore += lower; *ColCopy++ = *Mstore; }
      }
   }
}

