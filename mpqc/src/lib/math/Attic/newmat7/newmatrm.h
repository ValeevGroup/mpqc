//$$newmatrm.h                            rectangular matrix operations

// Copyright (C) 1991,2,3: R B Davies

#ifndef MATRIXRM_LIB
#define MATRIXRM_LIB 0


// operations on rectangular matrices

class RectMatrixCol;

class RectMatrixRowCol
// a class for accessing rows and columns of rectangular matrices
{
protected:
   Real* store;                   // pointer to storage
   int n;                         // number of elements
   int spacing;                   // space between elements
   int shift;                     // space between cols or rows
   RectMatrixRowCol(Real* st, int nx, int sp, int sh);
   void Reset(Real* st, int nx, int sp, int sh);
public:
   Real operator*(const RectMatrixRowCol&) const;         // dot product
   void AddScaled(const RectMatrixRowCol&, Real);         // add scaled
   void Divide(const RectMatrixRowCol&, Real);            // scaling
   void Divide(Real);                                     // scaling
   void Negate();                                         // change sign
   void Zero();                                           // zero row col
   Real& operator[](int i); // element
   Real SumSquare() const;                                // sum of squares
   Real& First();                       // get first element
   void DownDiag();
   void UpDiag();
   friend void ComplexScale(RectMatrixCol&, RectMatrixCol&, Real, Real);
   friend void Rotate(RectMatrixCol&, RectMatrixCol&, Real, Real);
   FREE_CHECK(RectMatrixRowCol)
};

class RectMatrixRow : public RectMatrixRowCol
{
public:
   RectMatrixRow(const Matrix&, int, int, int);
   RectMatrixRow(const Matrix&, int);
   void Reset(const Matrix&, int, int, int);
   void Reset(const Matrix&, int);
   Real& operator[](int i);
   void Down();
   void Right();
   void Up();
   void Left();
   FREE_CHECK(RectMatrixRow)
};

class RectMatrixCol : public RectMatrixRowCol
{
public:
   RectMatrixCol(const Matrix&, int, int, int);
   RectMatrixCol(const Matrix&, int);
   void Reset(const Matrix&, int, int, int);
   void Reset(const Matrix&, int);
   void Down();
   void Right();
   void Up();
   void Left();
   friend void ComplexScale(RectMatrixCol&, RectMatrixCol&, Real, Real);
   friend void Rotate(RectMatrixCol&, RectMatrixCol&, Real, Real);
   FREE_CHECK(RectMatrixCol)
};

class RectMatrixDiag : public RectMatrixRowCol
{
public:
   RectMatrixDiag(const DiagonalMatrix& D);
   Real& operator[](int i);
   void DownDiag();
   void UpDiag();
   FREE_CHECK(RectMatrixDiag)
};


 Real square(Real x);
 Real sign(Real x, Real y);


#endif
