//$$ newmatrc.h              definition file for row/column classes

// Copyright (C) 1991,2,3: R B Davies

#ifndef MATRIXRC_LIB
#define MATRIXRC_LIB 0

#include <math/newmat7/controlw.h>



/************** classes MatrixRowCol, MatrixRow, MatrixCol *****************/

// Used for accessing the rows and columns of matrices
// All matrix classes must provide routines for calculating matrix rows and
// columns. Assume rows can be found very efficiently.

enum LSF { LoadOnEntry=1,StoreOnExit=2,IsACopy=4,DirectPart=8,StoreHere=16 };


class LoadAndStoreFlag : public ControlWord
{
public:
   LoadAndStoreFlag();
   LoadAndStoreFlag(int i);
   LoadAndStoreFlag(LSF lsf);
   LoadAndStoreFlag(const ControlWord& cwx);
   FREE_CHECK(LoadAndStoreFlag)
};

class MatrixRowCol
// the row or column of a matrix
{
public:                                        // these are public to avoid
					       // numerous friend statements
   int length;                                 // row or column length
   int skip;                                   // initial number of zeros
   int storage;                                // number of stored elements
   int rowcol;                                 // row or column number
   GeneralMatrix* gm;                          // pointer to parent matrix
   Real* store;                                // pointer to local storage
					       //    less skip
   LoadAndStoreFlag cw;                        // Load? Store? Is a Copy?
   void IncrMat();
					       // used by NextRow
   void IncrDiag();
   void IncrUT();
   void IncrLT();

public:
   void Add(const MatrixRowCol&);              // add a row/col
   void AddScaled(const MatrixRowCol&, Real);  // add a multiple of a row/col
   void Add(const MatrixRowCol&, const MatrixRowCol&);
					       // add two rows/cols
   void Add(const MatrixRowCol&, Real);        // add a row/col
   void Sub(const MatrixRowCol&);              // subtract a row/col
   void Sub(const MatrixRowCol&, const MatrixRowCol&);
					       // sub a row/col from another
   void RevSub(const MatrixRowCol&);           // subtract from a row/col
   void Copy(const MatrixRowCol&);             // copy a row/col
   void CopyCheck(const MatrixRowCol&);        // ... check for data loss
   void Copy(const Real*&);                    // copy from an array
   void Copy(Real);                            // copy from constant
   Real SumAbsoluteValue();                    // sum of absolute values
   void Inject(const MatrixRowCol&);           // copy stored els of a row/col
   void Negate(const MatrixRowCol&);           // change sign of a row/col
   void Multiply(const MatrixRowCol&, Real);   // scale a row/col
   friend Real DotProd(const MatrixRowCol&, const MatrixRowCol&);
					       // sum of pairwise product
   Real* operator()();   // pointer to first element
   Real* Store();
   int Skip();                 // number of elements skipped
   int Storage();           // number of elements stored
   void Skip(int i);
   void Storage(int i);
   void SubRowCol(MatrixRowCol&, int, int) const;
					       // get part of a row or column
   MatrixRowCol();                           // to stop warning messages
   ~MatrixRowCol();
   FREE_CHECK(MatrixRowCol)
};

class MatrixRow : public MatrixRowCol
{
public:
   // bodies for these are inline at the end of this .h file
   MatrixRow(GeneralMatrix*, LoadAndStoreFlag, int=0);
                                               // extract a row
   ~MatrixRow();
   void Next();                                // get next row
   FREE_CHECK(MatrixRow)
};

class MatrixCol : public MatrixRowCol
{
public:
   // bodies for these are inline at the end of this .h file
   MatrixCol(GeneralMatrix*, LoadAndStoreFlag, int=0);
                                               // extract a col
   MatrixCol(GeneralMatrix*, Real*, LoadAndStoreFlag, int=0);
                                               // store/retrieve a col
   ~MatrixCol();
   void Next();                                // get next row
   FREE_CHECK(MatrixCol)
};


/**************************** inline bodies ****************************/


#endif
