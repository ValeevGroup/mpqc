//$$ newmat.h           definition file for new version of matrix package

// Copyright (C) 1991,2,3: R B Davies

#ifndef MATRIX_LIB
#define MATRIX_LIB 0

#include <math.h>
#ifdef NO_LONG_NAMES
#define UpperTriangularMatrix UTMatrix
#define LowerTriangularMatrix LTMatrix
#define SymmetricMatrix SMatrix
#define DiagonalMatrix DMatrix
#endif

#include <math/newmat7/include.h>

#ifndef TEMPS_DESTROYED_QUICKLY
#define ReturnMatrix ReturnMatrixX
#else
#define ReturnMatrix ReturnMatrixX&
#endif

#include <math/newmat7/boolean.h>
#include <math/newmat7/except.h>

/**************************** general utilities ****************************/

class GeneralMatrix;
void MatrixErrorNoSpace(void*);                 // no space handler

class LogAndSign
// Return from LogDeterminant function
//    - value of the log plus the sign (+, - or 0)
{
   Real log_value;
   int sign;
public:
   LogAndSign();
   LogAndSign(Real);
   void operator*=(Real);
   void ChangeSign();
   Real LogValue() const;
   int Sign() const;
   Real Value() const;
   FREE_CHECK(LogAndSign)
};

// the following class is for counting the number of times a piece of code
// is executed. It is used for locating any code not executed by test
// routines. Use turbo GREP locate all places this code is called and
// check which ones are not accessed.
// Somewhat implementation dependent as it relies on "cout" still being
// present when ExeCounter objects are destructed.

class ExeCounter
{
   int line;                                    // code line number
   int fileid;                                  // file identifier
   long nexe;                                   // number of executions
   static int nreports;                         // number of reports
public:
   ExeCounter(int,int);
   void operator++();
   ~ExeCounter();                               // prints out reports
};


/************** Class to show whether to check for loss of data ************/

class MatrixConversionCheck : public Janitor
{
   static Boolean DoCheck;
public:
   MatrixConversionCheck();          // turn check on
   ~MatrixConversionCheck();        // turn check off
   void CleanUp();
   static Boolean IsOn();
   static void DataLoss();
};

/**************************** class MatrixType *****************************/

// Is used for finding the type of a matrix resulting from the binary operations
// +, -, * and identifying what conversions are permissible.
// This class must be updated when new matrix types are added.

class GeneralMatrix;                            // defined later
class BaseMatrix;                               // defined later
class MatrixInput;                              // defined later

class MatrixType
{
public:
   enum Attribute {  Valid     = 1,
                     Symmetric = 2,
                     Band      = 4,
                     Upper     = 8,
                     Lower     = 16,
                     LUDeco    = 32 };

   enum            { US = 0,
                     UT = Valid + Upper,
                     LT = Valid + Lower,
                     Rt = Valid,
                     Sm = Valid + Symmetric,
                     Dg = Valid + Band + Lower + Upper + Symmetric,
		     RV = Valid,     //   don't separate out
		     CV = Valid,     //   vectors
		     BM = Valid + Band,
		     UB = Valid + Band + Upper,
		     LB = Valid + Band + Lower,
		     SB = Valid + Band + Symmetric,
		     Ct = Valid + LUDeco,
		     BC = Valid + Band + LUDeco,
                   };


   static nTypes();               // number of different types
					       // exclude Ct, US, BC
public:
   int attribute;
public:
   MatrixType ();
   MatrixType (int i);
   int operator+() const;
   MatrixType operator+(MatrixType mt) const;
   MatrixType operator*(const MatrixType&) const;
   Boolean operator>=(MatrixType mt) const;
   Boolean operator==(MatrixType t) const;
   Boolean operator!=(MatrixType t) const;
   Boolean operator!() const;
   MatrixType i() const;                        // type of inverse
   MatrixType t() const;                       // type of transpose
   MatrixType AddEqualEl() const;               // Add constant to matrix
   MatrixType MultRHS() const;                 // type for rhs of multiply
   MatrixType sub() const;                      // type of submatrix
   MatrixType ssub() const;                     // type of sym submatrix
   GeneralMatrix* New() const;                 // new matrix of given type
   GeneralMatrix* New(int,int,BaseMatrix*) const;
                                               // new matrix of given type
   operator char*() const;                     // to print type
   friend Boolean Rectangular(MatrixType a, MatrixType b, MatrixType c);
//   friend Boolean Rectangular(MatrixType a, MatrixType b, MatrixType c)
//      { return ((a.attribute | b.attribute | c.attribute)
//          & ~MatrixType::Valid) == 0; }
   friend Boolean Compare(const MatrixType&, MatrixType&);
                                               // compare and check conv.
   ~MatrixType();
   FREE_CHECK(MatrixType)
};

void TestTypeAdd();                            // test +
void TestTypeMult();                           // test *
void TestTypeOrder();                          // test >=


/************************* class MatrixBandWidth ***********************/

class MatrixBandWidth
{
public:
   int lower;
   int upper;
   MatrixBandWidth(const int l, const int u);
   MatrixBandWidth(const int i);
   MatrixBandWidth operator+(const MatrixBandWidth&) const;
   MatrixBandWidth operator*(const MatrixBandWidth&) const;
   MatrixBandWidth t() const;
   Boolean operator==(const MatrixBandWidth& bw) const;
   FREE_CHECK(MatrixBandWidth)
};


/*********************** Array length specifier ************************/

// This class is introduced to avoid constructors such as
//   ColumnVector(int)
// being used for conversions

class ArrayLengthSpecifier
{
   int value;
public:
   int Value() const;
   ArrayLengthSpecifier(int l);
   FREE_CHECK(ArrayLengthSpecifier)
};

/*************************** Matrix routines ***************************/


class MatrixRowCol;                             // defined later
class MatrixRow;
class MatrixCol;

class GeneralMatrix;                            // defined later
class AddedMatrix;
class MultipliedMatrix;
class SubtractedMatrix;
class SolvedMatrix;
class ShiftedMatrix;
class ScaledMatrix;
class TransposedMatrix;
class NegatedMatrix;
class InvertedMatrix;
class RowedMatrix;
class ColedMatrix;
class DiagedMatrix;
class MatedMatrix;
class GetSubMatrix;
class ConstMatrix;
class ReturnMatrixX;
class Matrix;
class nricMatrix;
class RowVector;
class ColumnVector;
class SymmetricMatrix;
class UpperTriangularMatrix;
class LowerTriangularMatrix;
class DiagonalMatrix;
class CroutMatrix;
class BandMatrix;
class LowerBandMatrix;
class UpperBandMatrix;
class SymmetricBandMatrix;
class LinearEquationSolver;



static MatrixType MatrixTypeUnSp(MatrixType::US);
						// AT&T needs this

class BaseMatrix : public Janitor               // base of all matrix classes
{
protected:
   virtual int search(const BaseMatrix*) const = 0;
						// count number of times matrix
						// is referred to
public:
   ~BaseMatrix();
#ifndef __GNUG__
   virtual GeneralMatrix* Evaluate(MatrixType mt=MatrixTypeUnSp) = 0;
						// evaluate temporary
#else
   virtual GeneralMatrix* Evaluate(MatrixType mt) = 0;
   GeneralMatrix* Evaluate();
#endif
#ifndef TEMPS_DESTROYED_QUICKLY
   AddedMatrix operator+(const BaseMatrix&) const;    // results of operations
   MultipliedMatrix operator*(const BaseMatrix&) const;
   SubtractedMatrix operator-(const BaseMatrix&) const;
   ShiftedMatrix operator+(Real) const;
   ScaledMatrix operator*(Real) const;
   ScaledMatrix operator/(Real) const;
   ShiftedMatrix operator-(Real) const;
   TransposedMatrix t() const;
//   TransposedMatrix t;
   NegatedMatrix operator-() const;                   // change sign of elements
   InvertedMatrix i() const;
//   InvertedMatrix i;
   RowedMatrix AsRow() const;
   ColedMatrix AsColumn() const;
   DiagedMatrix AsDiagonal() const;
   MatedMatrix AsMatrix(int,int) const;
   GetSubMatrix SubMatrix(int,int,int,int) const;
   GetSubMatrix SymSubMatrix(int,int) const;
   GetSubMatrix Row(int) const;
   GetSubMatrix Rows(int,int) const;
   GetSubMatrix Column(int) const;
   GetSubMatrix Columns(int,int) const;
#else
   AddedMatrix& operator+(const BaseMatrix&) const;    // results of operations
   MultipliedMatrix& operator*(const BaseMatrix&) const;
   SubtractedMatrix& operator-(const BaseMatrix&) const;
   ShiftedMatrix& operator+(Real) const;
   ScaledMatrix& operator*(Real) const;
   ScaledMatrix& operator/(Real) const;
   ShiftedMatrix& operator-(Real) const;
   TransposedMatrix& t() const;
//   TransposedMatrix& t;
   NegatedMatrix& operator-() const;                   // change sign of elements
   InvertedMatrix& i() const;
//   InvertedMatrix& i;
   RowedMatrix& AsRow() const;
   ColedMatrix& AsColumn() const;
   DiagedMatrix& AsDiagonal() const;
   MatedMatrix& AsMatrix(int,int) const;
   GetSubMatrix& SubMatrix(int,int,int,int) const;
   GetSubMatrix& SymSubMatrix(int,int) const;
   GetSubMatrix& Row(int) const;
   GetSubMatrix& Rows(int,int) const;
   GetSubMatrix& Column(int) const;
   GetSubMatrix& Columns(int,int) const;
#endif
   Real AsScalar() const;                      // conversion of 1 x 1 matrix
   virtual LogAndSign LogDeterminant() const;
   virtual Real SumSquare() const;
   virtual Real SumAbsoluteValue() const;
   virtual Real MaximumAbsoluteValue() const;
   virtual Real Trace() const;
   Real Norm1() const;
   Real NormInfinity() const;
   virtual MatrixBandWidth BandWidth() const;  // bandwidths of band matrix
   virtual void CleanUp();                     // to clear store
//protected:
//   BaseMatrix() : t(this), i(this) {}

   friend class GeneralMatrix;
   friend class Matrix;
   friend class nricMatrix;
   friend class RowVector;
   friend class ColumnVector;
   friend class SymmetricMatrix;
   friend class UpperTriangularMatrix;
   friend class LowerTriangularMatrix;
   friend class DiagonalMatrix;
   friend class CroutMatrix;
   friend class BandMatrix;
   friend class LowerBandMatrix;
   friend class UpperBandMatrix;
   friend class SymmetricBandMatrix;
   friend class AddedMatrix;
   friend class MultipliedMatrix;
   friend class SubtractedMatrix;
   friend class SolvedMatrix;
   friend class ShiftedMatrix;
   friend class ScaledMatrix;
   friend class TransposedMatrix;
   friend class NegatedMatrix;
   friend class InvertedMatrix;
   friend class RowedMatrix;
   friend class ColedMatrix;
   friend class DiagedMatrix;
   friend class MatedMatrix;
   friend class GetSubMatrix;
   friend class ConstMatrix;
   friend class ReturnMatrixX;
   friend class LinearEquationSolver;
   NEW_DELETE(BaseMatrix)
};


/******************************* working classes **************************/

class GeneralMatrix : public BaseMatrix         // declarable matrix types
{
   virtual GeneralMatrix* Image() const;        // copy of matrix
protected:
   int tag;                                     // shows whether can reuse
   int nrows, ncols;                            // dimensions
   int storage;                                 // total store required
   Real* store;                                 // point to store (0=not set)
   GeneralMatrix();                             // initialise with no store
   GeneralMatrix(ArrayLengthSpecifier);         // constructor getting store
   void Add(GeneralMatrix*, Real);              // sum of GM and Real
   void Add(Real);                              // add Real to this
   void Multiply(GeneralMatrix*, Real);         // product of GM and Real
   void Multiply(Real);                         // multiply this by Real
   void Negate(GeneralMatrix*);                 // change sign
   void Negate();                               // change sign
   void operator=(Real);                        // set matrix to constant
   Real* GetStore();                            // get store or copy
   GeneralMatrix* BorrowStore(GeneralMatrix*, MatrixType);
                                                // temporarily access store
   void GetMatrix(const GeneralMatrix*);        // used by = and initialise
   void Eq(const BaseMatrix&, MatrixType);      // used by =
   int search(const BaseMatrix*) const;
   virtual GeneralMatrix* Transpose(TransposedMatrix*, MatrixType);
   void CheckConversion(const BaseMatrix&);     // check conversion OK
   void ReDimension(int, int, int);             // change dimensions
public:
   GeneralMatrix* Evaluate(MatrixType);
   virtual MatrixType Type() const = 0;       // type of a matrix
   int Nrows() const;          // get dimensions
   int Ncols() const;
   int Storage() const;
   Real* Store() const;
   virtual ~GeneralMatrix();                    // delete store if set
   void tDelete();                              // delete if tag permits
   Boolean reuse();                                // TRUE if tag allows reuse
   void Protect();                   // can't delete or reuse
   int Tag() const;
   Boolean IsZero() const;                         // test matrix has all zeros
   void Release();                    // del store after next use
   void Release(int t);               // del store after t accesses
   void ReleaseAndDelete();           // delete matrix after use
   void operator<<(const Real*);                // assignment from an array
   void operator<<(const BaseMatrix& X);
                                                // = without checking type
   void Inject(const GeneralMatrix&);           // copy stored els only
   virtual GeneralMatrix* MakeSolver();         // for solving
   virtual void Solver(MatrixRowCol&, const MatrixRowCol&);
   virtual void GetRow(MatrixRowCol&) = 0;      // Get matrix row
   virtual void RestoreRow(MatrixRowCol&);    // Restore matrix row
   virtual void NextRow(MatrixRowCol&);         // Go to next row
   virtual void GetCol(MatrixRowCol&) = 0;      // Get matrix col
   virtual void RestoreCol(MatrixRowCol&);    // Restore matrix col
   virtual void NextCol(MatrixRowCol&);         // Go to next col
   Real SumSquare() const;
   Real SumAbsoluteValue() const;
   Real MaximumAbsoluteValue() const;
   LogAndSign LogDeterminant() const;
#ifndef TEMPS_DESTROYED_QUICKLY
   ConstMatrix c() const;                       // to access constant matrices
#else
   ConstMatrix& c() const;                      // to access constant matrices
#endif
   void CheckStore() const;                     // check store is non-zero
   virtual void SetParameters(const GeneralMatrix*);
                                                // set parameters in GetMatrix
#ifndef TEMPS_DESTROYED_QUICKLY
   operator ReturnMatrixX() const;              // for building a ReturnMatrix
#else
   operator ReturnMatrixX&() const;             // for building a ReturnMatrix
#endif
   MatrixInput operator<<(Real);                // for loading a list
   void CleanUp();                              // to clear store

   friend class Matrix;
   friend class nricMatrix;
   friend class SymmetricMatrix;
   friend class UpperTriangularMatrix;
   friend class LowerTriangularMatrix;
   friend class DiagonalMatrix;
   friend class CroutMatrix;
   friend class RowVector;
   friend class ColumnVector;
   friend class BandMatrix;
   friend class LowerBandMatrix;
   friend class UpperBandMatrix;
   friend class SymmetricBandMatrix;
   friend class BaseMatrix;
   friend class AddedMatrix;
   friend class MultipliedMatrix;
   friend class SubtractedMatrix;
   friend class SolvedMatrix;
   friend class ShiftedMatrix;
   friend class ScaledMatrix;
   friend class TransposedMatrix;
   friend class NegatedMatrix;
   friend class InvertedMatrix;
   friend class RowedMatrix;
   friend class ColedMatrix;
   friend class DiagedMatrix;
   friend class MatedMatrix;
   friend class GetSubMatrix;
   friend class ConstMatrix;
   friend class ReturnMatrixX;
   friend class LinearEquationSolver;
   NEW_DELETE(GeneralMatrix)
};

class Matrix : public GeneralMatrix             // usual rectangular matrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~Matrix();
   Matrix();
   Matrix(int, int);                            // standard declaration
   Matrix(const BaseMatrix&);                   // evaluate BaseMatrix
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const Matrix& m);
   MatrixType Type() const;
   Real& operator()(int, int);                  // access element
   Real& element(int, int);                     // access element
   const Real& element(int, int) const;                     // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   Matrix(const Matrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
#endif
   GeneralMatrix* MakeSolver();
   Real Trace() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void RestoreCol(MatrixRowCol&);
   void NextRow(MatrixRowCol&);
   void NextCol(MatrixRowCol&);
   void ReDimension(int,int);                   // change dimensions
   NEW_DELETE(Matrix)
};

class nricMatrix : public Matrix                // for use with Numerical
                                                // Recipes in C
{
   GeneralMatrix* Image() const;                // copy of matrix
   Real** row_pointer;                          // points to rows
   void MakeRowPointer();                       // build rowpointer
   void DeleteRowPointer();
public:
   nricMatrix();
   nricMatrix(int m, int n);                     // standard declaration
   nricMatrix(const BaseMatrix& bm);             // evaluate BaseMatrix
   void operator=(const BaseMatrix& bm);
   void operator=(Real f);
   void operator=(const nricMatrix& m);
   void operator<<(const BaseMatrix& X);
   nricMatrix(const nricMatrix& gm);
   void ReDimension(int m, int n);               // change dimensions
   ~nricMatrix();
#ifndef __ZTC__
   Real** nric() const;
#endif
   void CleanUp();                                // to clear store
   NEW_DELETE(nricMatrix)
};

class SymmetricMatrix : public GeneralMatrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~SymmetricMatrix();
   SymmetricMatrix();
   SymmetricMatrix(ArrayLengthSpecifier);
   SymmetricMatrix(const BaseMatrix&);
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const SymmetricMatrix& m);
   Real& operator()(int, int);                  // access element
   Real& element(int, int);                     // access element
   const Real& element(int, int) const;                     // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   MatrixType Type() const;
   SymmetricMatrix(const SymmetricMatrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
#endif
   Real SumSquare() const;
   Real SumAbsoluteValue() const;
   Real Trace() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   GeneralMatrix* Transpose(TransposedMatrix*, MatrixType);
   void ReDimension(int);                       // change dimensions
   NEW_DELETE(SymmetricMatrix)
};

class UpperTriangularMatrix : public GeneralMatrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~UpperTriangularMatrix();
   UpperTriangularMatrix();
   UpperTriangularMatrix(ArrayLengthSpecifier);
   void operator=(const BaseMatrix&);
   void operator=(const UpperTriangularMatrix& m);
   UpperTriangularMatrix(const BaseMatrix&);
   UpperTriangularMatrix(const UpperTriangularMatrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
#endif
   void operator=(Real f);
   Real& operator()(int, int);                  // access element
   Real& element(int, int);                     // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   MatrixType Type() const;
   GeneralMatrix* MakeSolver(); // for solving
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   LogAndSign LogDeterminant() const;
   Real Trace() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void RestoreCol(MatrixRowCol&);
   void NextRow(MatrixRowCol&);
   void ReDimension(int);                       // change dimensions
   NEW_DELETE(UpperTriangularMatrix)
};

class LowerTriangularMatrix : public GeneralMatrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~LowerTriangularMatrix();
   LowerTriangularMatrix();
   LowerTriangularMatrix(ArrayLengthSpecifier);
   LowerTriangularMatrix(const LowerTriangularMatrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
#endif
   LowerTriangularMatrix(const BaseMatrix& M);
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const LowerTriangularMatrix& m);
   Real& operator()(int, int);                  // access element
   Real& element(int, int);                     // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   MatrixType Type() const;
   GeneralMatrix* MakeSolver(); // for solving
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   LogAndSign LogDeterminant() const;
   Real Trace() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void RestoreCol(MatrixRowCol&);
   void NextRow(MatrixRowCol&);
   void ReDimension(int);                       // change dimensions
   NEW_DELETE(LowerTriangularMatrix)
};

class DiagonalMatrix : public GeneralMatrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~DiagonalMatrix();
   DiagonalMatrix();
   DiagonalMatrix(ArrayLengthSpecifier);
   DiagonalMatrix(const BaseMatrix&);
   DiagonalMatrix(const DiagonalMatrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
   Real operator()(int) const;
#endif
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const DiagonalMatrix& m);
   Real& operator()(int, int);                  // access element
   Real& operator()(int);                       // access element
   Real& element(int, int);                     // access element
   Real& element(int);                          // access element
   const Real& element(int, int) const;                    // access element
   const Real& element(int) const;                          // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real& operator[](int m);
   const Real& operator[](int m) const;
#endif
   MatrixType Type() const;

   LogAndSign LogDeterminant() const;
   Real Trace() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void NextRow(MatrixRowCol&);
   void NextCol(MatrixRowCol&);
   GeneralMatrix* MakeSolver(); // for solving
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   GeneralMatrix* Transpose(TransposedMatrix*, MatrixType);
   void ReDimension(int);                       // change dimensions
#ifndef __ZTC__
   Real* nric() const;
#endif
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(DiagonalMatrix)
};

class RowVector : public Matrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~RowVector();
   RowVector();
   RowVector(ArrayLengthSpecifier n);
   RowVector(const BaseMatrix&);
   RowVector(const RowVector& gm);
#ifndef __ZTC__
   Real operator()(int) const;                  // access element
#endif
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const RowVector& m);
   Real& operator()(int);                       // access element
   Real& element(int);                          // access element
   const Real& element(int) const;              // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real& operator[](int m);
   const Real& operator[](int m) const;
#endif
   MatrixType Type() const;
   void GetCol(MatrixRowCol&);
   void NextCol(MatrixRowCol&);
   void RestoreCol(MatrixRowCol&);
   GeneralMatrix* Transpose(TransposedMatrix*, MatrixType);
   void ReDimension(int);                       // change dimensions
#ifndef __ZTC__
   Real* nric() const;
#endif
   void CleanUp();                                // to clear store
   NEW_DELETE(RowVector)
};

class ColumnVector : public Matrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~ColumnVector();
   ColumnVector();
   ColumnVector(ArrayLengthSpecifier n);
   ColumnVector(const BaseMatrix&);
   ColumnVector(const ColumnVector& gm);
#ifndef __ZTC__
   Real operator()(int) const;                  // access element
#endif
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const ColumnVector& m);
   Real& operator()(int);                       // access element
   Real& element(int);                          // access element
   const Real& element(int) const;              // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real& operator[](int m);
   const Real& operator[](int m) const;
#endif
   MatrixType Type() const;
   GeneralMatrix* Transpose(TransposedMatrix*, MatrixType);
   void ReDimension(int);                       // change dimensions
#ifndef __ZTC__
   Real* nric() const;
#endif
   void CleanUp();                                // to clear store
   NEW_DELETE(ColumnVector)
//JCM
   Real Dot(ColumnVector& y);
   Real Norm2(ColumnVector& y);
};

class CroutMatrix : public GeneralMatrix        // for LU decomposition
{
   int* indx;
   Boolean d;
   Boolean sing;
   void ludcmp();
public:
   CroutMatrix(const BaseMatrix&);
   MatrixType Type() const;
   void lubksb(Real*, int=0);
   ~CroutMatrix();
   GeneralMatrix* MakeSolver(); // for solving
   LogAndSign LogDeterminant() const;
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void operator=(const BaseMatrix&);
   void operator=(const CroutMatrix& m);
   void CleanUp();                                // to clear store
   NEW_DELETE(CroutMatrix)
};

/******************************* band matrices ***************************/

class BandMatrix : public GeneralMatrix         // band matrix
{
   GeneralMatrix* Image() const;                // copy of matrix
protected:
   void CornerClear() const;                    // set unused elements to zero
public:
   ~BandMatrix();
   int lower, upper;                            // band widths
   BandMatrix();
   BandMatrix(int n,int lb,int ub);
                                                // standard declaration
   BandMatrix(const BaseMatrix&);               // evaluate BaseMatrix
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const BandMatrix& m);
   MatrixType Type() const;
   Real& operator()(int, int);                  // access element
   Real& element(int, int);                     // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   BandMatrix(const BandMatrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
#endif
   LogAndSign LogDeterminant() const;
   GeneralMatrix* MakeSolver();
   Real Trace() const;
   Real SumSquare() const;
   Real SumAbsoluteValue() const;
   Real MaximumAbsoluteValue() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void RestoreCol(MatrixRowCol&);
   void NextRow(MatrixRowCol&);
   void ReDimension(int, int, int);             // change dimensions
   MatrixBandWidth BandWidth() const;
   void SetParameters(const GeneralMatrix*);
   MatrixInput operator<<(Real);                // will give error
   void operator<<(const Real* r);              // will give error
      // the next is included because Zortech and Borland
      // can't find the copy in GeneralMatrix
   void operator<<(const BaseMatrix& X);
   NEW_DELETE(BandMatrix)
};

class UpperBandMatrix : public BandMatrix       // upper band matrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~UpperBandMatrix();
   UpperBandMatrix();
   UpperBandMatrix(int n, int ubw);              // standard declaration
   UpperBandMatrix(const BaseMatrix&);          // evaluate BaseMatrix
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const UpperBandMatrix& m);
   MatrixType Type() const;
   UpperBandMatrix(const UpperBandMatrix& gm);
   GeneralMatrix* MakeSolver();
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   LogAndSign LogDeterminant() const;
   void ReDimension(int n,int ubw);              // change dimensions
   Real& operator()(int, int);
   Real& element(int, int);
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   NEW_DELETE(UpperBandMatrix)
};

class LowerBandMatrix : public BandMatrix       // upper band matrix
{
   GeneralMatrix* Image() const;                // copy of matrix
public:
   ~LowerBandMatrix();
   LowerBandMatrix();
   LowerBandMatrix(int n, int lbw);              // standard declaration
   LowerBandMatrix(const BaseMatrix&);          // evaluate BaseMatrix
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const LowerBandMatrix& m);
   MatrixType Type() const;
   LowerBandMatrix(const LowerBandMatrix& gm);
   GeneralMatrix* MakeSolver();
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   LogAndSign LogDeterminant() const;
   void ReDimension(int n,int lbw);             // change dimensions
   Real& operator()(int, int);
   Real& element(int, int);
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   NEW_DELETE(LowerBandMatrix)
};

class SymmetricBandMatrix : public GeneralMatrix
{
   GeneralMatrix* Image() const;                // copy of matrix
   void CornerClear() const;                    // set unused elements to zero
public:
   ~SymmetricBandMatrix();
   int lower;                                   // lower band width
   SymmetricBandMatrix();
   SymmetricBandMatrix(int n, int lb);
   SymmetricBandMatrix(const BaseMatrix&);
   void operator=(const BaseMatrix&);
   void operator=(Real f);
   void operator=(const SymmetricBandMatrix& m);
   Real& operator()(int, int);                  // access element
   Real& element(int, int);                     // access element
#ifdef SETUP_C_SUBSCRIPTS
   Real* operator[](int m);
   const Real* operator[](int m) const;
#endif
   MatrixType Type() const;
   SymmetricBandMatrix(const SymmetricBandMatrix& gm);
#ifndef __ZTC__
   Real operator()(int, int) const;             // access element
#endif
   GeneralMatrix* MakeSolver();
   Real SumSquare() const;
   Real SumAbsoluteValue() const;
   Real MaximumAbsoluteValue() const;
   Real Trace() const;
   LogAndSign LogDeterminant() const;
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   GeneralMatrix* Transpose(TransposedMatrix*, MatrixType);
   void ReDimension(int,int);                       // change dimensions
   MatrixBandWidth BandWidth() const;
   void SetParameters(const GeneralMatrix*);
   NEW_DELETE(SymmetricBandMatrix)
};

class BandLUMatrix : public GeneralMatrix
// for LU decomposition of band matrix
{
   int* indx;
   Boolean d;
   Boolean sing;                                   // TRUE if singular
   Real* store2;
   int storage2;
   void ludcmp();
   int m1,m2;                                   // lower and upper
public:
   BandLUMatrix(const BaseMatrix&);
   MatrixType Type() const;
   void lubksb(Real*, int=0);
   ~BandLUMatrix();
   GeneralMatrix* MakeSolver(); // for solving
   LogAndSign LogDeterminant() const;
   void Solver(MatrixRowCol&, const MatrixRowCol&);
   void GetRow(MatrixRowCol&);
   void GetCol(MatrixRowCol&);
   void operator=(const BaseMatrix&);
   void operator=(const BandLUMatrix& m);
   NEW_DELETE(BandLUMatrix)
   void CleanUp();                                // to clear store
};


/***************************** temporary classes *************************/

class MultipliedMatrix : public BaseMatrix
{
protected:
   union { const BaseMatrix* bm1; GeneralMatrix* gm1; };
						  // pointers to summands
   union { const BaseMatrix* bm2; GeneralMatrix* gm2; };
   MultipliedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x);
   int search(const BaseMatrix*) const;
   friend class BaseMatrix;
public:
   ~MultipliedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(MultipliedMatrix)
};

class AddedMatrix : public MultipliedMatrix
{
protected:
   AddedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x);

   friend class BaseMatrix;
public:
   ~AddedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
#ifdef __GNUG__
   void SelectVersion(MatrixType, int&, int&) const;
#else
   void SelectVersion(MatrixType, Boolean&, Boolean&) const;
#endif
   NEW_DELETE(AddedMatrix)
};

class SolvedMatrix : public MultipliedMatrix
{
   SolvedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x);
   friend class BaseMatrix;
   friend class InvertedMatrix;                        // for operator*
public:
   ~SolvedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(SolvedMatrix)
};

class SubtractedMatrix : public AddedMatrix
{
   SubtractedMatrix(const BaseMatrix* bm1x, const BaseMatrix* bm2x);
   friend class BaseMatrix;
public:
   ~SubtractedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   NEW_DELETE(SubtractedMatrix)
};

class ShiftedMatrix : public BaseMatrix
{
protected:
   Real f;
   union { const BaseMatrix* bm; GeneralMatrix* gm; };
   ShiftedMatrix(const BaseMatrix* bmx, Real fx);
   int search(const BaseMatrix*) const;
   friend class BaseMatrix;
public:
   ~ShiftedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   NEW_DELETE(ShiftedMatrix)
};

class ScaledMatrix : public ShiftedMatrix
{
   ScaledMatrix(const BaseMatrix* bmx, Real fx);
   friend class BaseMatrix;
public:
   ~ScaledMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(ScaledMatrix)
};

class NegatedMatrix : public BaseMatrix
{
protected:
   union { const BaseMatrix* bm; GeneralMatrix* gm; };
   NegatedMatrix(const BaseMatrix* bmx);
   int search(const BaseMatrix*) const;
private:
   friend class BaseMatrix;
public:
   ~NegatedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(NegatedMatrix)
};

class TransposedMatrix : public NegatedMatrix
{
   TransposedMatrix(const BaseMatrix* bmx);
   friend class BaseMatrix;
public:
   ~TransposedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(TransposedMatrix)
};

class InvertedMatrix : public NegatedMatrix
{
   InvertedMatrix(const BaseMatrix* bmx);
public:
   ~InvertedMatrix();
#ifndef TEMPS_DESTROYED_QUICKLY
   SolvedMatrix operator*(const BaseMatrix&) const;  // inverse(A) * B
#else
   SolvedMatrix& operator*(const BaseMatrix&);       // inverse(A) * B
#endif
   friend class BaseMatrix;
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(InvertedMatrix)
};

class RowedMatrix : public NegatedMatrix
{
   RowedMatrix(const BaseMatrix* bmx);
   friend class BaseMatrix;
public:
   ~RowedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(RowedMatrix)
};

class ColedMatrix : public NegatedMatrix
{
   ColedMatrix(const BaseMatrix* bmx);
   friend class BaseMatrix;
public:
   ~ColedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(ColedMatrix)
};

class DiagedMatrix : public NegatedMatrix
{
   DiagedMatrix(const BaseMatrix* bmx);
   friend class BaseMatrix;
public:
   ~DiagedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(DiagedMatrix)
};

class MatedMatrix : public NegatedMatrix
{
   int nr, nc;
   MatedMatrix(const BaseMatrix* bmx, int nrx, int ncx);
   friend class BaseMatrix;
public:
   ~MatedMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(MatedMatrix)
};

class ConstMatrix : public BaseMatrix
{
   const GeneralMatrix* cgm;
   int search(const BaseMatrix*) const;
   ConstMatrix(const GeneralMatrix* cgmx);
   friend class BaseMatrix;
   friend class GeneralMatrix;
public:
   ~ConstMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(ConstMatrix)
};

class ReturnMatrixX : public BaseMatrix    // for matrix return
{
   GeneralMatrix* gm;
   int search(const BaseMatrix*) const;
public:
   ~ReturnMatrixX();
   GeneralMatrix* Evaluate(MatrixType);
   friend class BaseMatrix;
#ifdef TEMPS_DESTROYED_QUICKLY
   ReturnMatrixX(const ReturnMatrixX& tm);
#else
   ReturnMatrixX(const ReturnMatrixX& tm);
#endif
   ReturnMatrixX(const GeneralMatrix* gmx);
//   ReturnMatrixX(GeneralMatrix&);
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(ReturnMatrixX)
};


/**************************** submatrices ******************************/

class GetSubMatrix : public NegatedMatrix
{
   int row_skip;
   int row_number;
   int col_skip;
   int col_number;
   Boolean IsSym;

   GetSubMatrix
      (const BaseMatrix* bmx, int rs, int rn, int cs, int cn, Boolean is);
   GetSubMatrix(const GetSubMatrix& g);
   void SetUpLHS();
   friend class BaseMatrix;
public:
   ~GetSubMatrix();
   GeneralMatrix* Evaluate(MatrixType);
   void operator=(const BaseMatrix&);
   void operator=(const GetSubMatrix& m);
   void operator<<(const BaseMatrix&);
   void operator<<(const Real*);                // copy from array
   void operator=(Real);                        // copy from constant
   void Inject(const GeneralMatrix&);           // copy stored els only
   MatrixBandWidth BandWidth() const;
   NEW_DELETE(GetSubMatrix)
};

/**************************** matrix input *******************************/

class MatrixInput          // for reading a list of values into a matrix
                           // the difficult part is detecting a mismatch
                           // in the number of elements
{
   static int n;           // number values still to be read
   static Real* r;         // pointer to last thing read
   static int depth;       // number of objects of this class in existence
public:
   MatrixInput();
   MatrixInput(const MatrixInput&);
   ~MatrixInput();
   MatrixInput operator<<(Real);
                           // could return a reference if we don't have
                           // TEMPS_DESTROYED_QUICKLY
   friend GeneralMatrix;                             
};

/***************************** exceptions ********************************/

class MatrixDetails
{
   MatrixType type;
   int nrows;
   int ncols;
   int ubw;
   int lbw;
public:
   MatrixDetails(const GeneralMatrix& A);
   void PrintOut();
   ~MatrixDetails();
};
   


class SpaceException : public Exception
{
public:
   static long st_type();
   long type() const;
   static int action;
   SpaceException();
   static void SetAction(int a);
};


class MatrixException : public Exception
{
public:
   static long st_type();
   long type() const;
   MatrixException(int);
   MatrixException(int, const GeneralMatrix&);
   MatrixException(int, const GeneralMatrix&, const GeneralMatrix&);
};

class DataException : public MatrixException
{
public:
   static long st_type();
   long type() const;
   static int action;
   DataException(const GeneralMatrix& A);
  static void SetAction(int a);
};

class SingularException : public DataException
{
public:
   static long st_type();
   long type() const;
   SingularException(const GeneralMatrix& A);
};

class NPDException : public DataException     // Not positive definite
{
public:
   static long st_type();
   long type() const;
   NPDException(const GeneralMatrix&);
};

class ConvergenceException : public MatrixException
{
public:
   static long st_type();
   long type() const;
   static int action;
   ConvergenceException(const GeneralMatrix& A);
   static void SetAction(int a);
};

class ProgramException : public MatrixException
{
public:
   static long st_type();
   long type() const;
   static int action;
   ProgramException(char* c);
   ProgramException(char* c, const GeneralMatrix&);
   ProgramException(char* c, const GeneralMatrix&, const GeneralMatrix&);
   ProgramException();
   ProgramException(const GeneralMatrix&);
   static void SetAction(int a);
};

class IndexException : public ProgramException
{
public:
   static long st_type();
   long type() const;
   IndexException(int i, const GeneralMatrix& A);
   IndexException(int i, int j, const GeneralMatrix& A);
   // next two are for access via element function
   IndexException(int i, const GeneralMatrix& A, Boolean);
   IndexException(int i, int j, const GeneralMatrix& A, Boolean);
};

class VectorException : public ProgramException    // can't convert to vector
{
public:
   static long st_type();
   long type() const;
   VectorException();
   VectorException(const GeneralMatrix& A);
};

class NotSquareException : public ProgramException
{
public:
   static long st_type();
   long type() const;
   NotSquareException(const GeneralMatrix& A);
};

class SubMatrixDimensionException : public ProgramException
{
public:
   static long st_type();
   long type() const;
   SubMatrixDimensionException();
};

class IncompatibleDimensionsException : public ProgramException
{
public:
   static long st_type();
   long type() const;
   IncompatibleDimensionsException();
};

class NotDefinedException : public ProgramException
{
public:
   static long st_type();
   long type() const;
   NotDefinedException(char* op, char* matrix);
};

class CannotBuildException : public ProgramException
{
public:
   static long st_type();
   long type() const;
   CannotBuildException(char* matrix);
};


class InternalException : public MatrixException
{
public:
   static long st_type();
   long type() const;
   static int action;
   InternalException(char* c);
   static void SetAction(int a);
};


/***************************** functions ***********************************/


 LogAndSign LogDeterminant(const BaseMatrix& B);
 Real SumSquare(const BaseMatrix& B);
 Real Trace(const BaseMatrix& B);
 Real SumAbsoluteValue(const BaseMatrix& B);
 Real MaximumAbsoluteValue(const BaseMatrix& B);
 Real Norm1(const BaseMatrix& B);
 Real Norm1(RowVector& RV);
 Real NormInfinity(const BaseMatrix& B);
 Real NormInfinity(ColumnVector& CV);

//JCM

 Real Norm2(const ColumnVector& CV);
 Real Dot(ColumnVector& CV1, ColumnVector& CV2);


//------------------------------------------------------------------------
// Print various quantities in newmat classes
//------------------------------------------------------------------------

void Print(const Matrix& X);
void Print(const UpperTriangularMatrix& X);
void Print(const DiagonalMatrix& X);
void Print(const SymmetricMatrix& X);
void Print(const LowerTriangularMatrix& X);

#endif
