//$$ newmatap.h           definition file for matrix package applications

// Copyright (C) 1991,2,3: R B Davies

#ifndef MATRIXAP_LIB
#define MATRIXAP_LIB 0

#include <math/newmat7/newmat.h>


/**************************** applications *****************************/


void HHDecompose(Matrix&, LowerTriangularMatrix&);

void HHDecompose(const Matrix&, Matrix&, Matrix&);

ReturnMatrix Cholesky(const SymmetricMatrix&);

ReturnMatrix Cholesky(const SymmetricBandMatrix&);

#ifndef __GNUG__
void SVD(const Matrix&, DiagonalMatrix&, Matrix&, Matrix&,
    Boolean=TRUE, Boolean=TRUE);
#else
void SVD(const Matrix&, DiagonalMatrix&, Matrix&, Matrix&,
    Boolean=(Boolean)TRUE, Boolean=(Boolean)TRUE);
#endif


void SVD(const Matrix&, DiagonalMatrix&);

#ifndef __GNUG__
 void SVD(const Matrix& A, DiagonalMatrix& D, Matrix& U,
          Boolean withU = TRUE);
#else
 void SVD(const Matrix& A, DiagonalMatrix& D, Matrix& U,
          Boolean withU = (Boolean)TRUE);
#endif

void Jacobi(const SymmetricMatrix&, DiagonalMatrix&);

void Jacobi(const SymmetricMatrix&, DiagonalMatrix&, SymmetricMatrix&);

void Jacobi(const SymmetricMatrix&, DiagonalMatrix&, Matrix&);

#ifndef __GNUG__
void Jacobi(const SymmetricMatrix&, DiagonalMatrix&, SymmetricMatrix&,
   Matrix&, Boolean=TRUE);
#else
void Jacobi(const SymmetricMatrix&, DiagonalMatrix&, SymmetricMatrix&,
   Matrix&, Boolean=(Boolean)TRUE);
#endif

void EigenValues(const SymmetricMatrix&, DiagonalMatrix&);

void EigenValues(const SymmetricMatrix&, DiagonalMatrix&, SymmetricMatrix&);

void EigenValues(const SymmetricMatrix&, DiagonalMatrix&, Matrix&);

class SymmetricEigenAnalysis
// not implemented yet
{
public:
   ~SymmetricEigenAnalysis();
   SymmetricEigenAnalysis(const SymmetricMatrix&);
private:
   DiagonalMatrix diag;
   DiagonalMatrix offdiag;
   SymmetricMatrix backtransform;
   FREE_CHECK(SymmetricEigenAnalysis)
};

void SortAscending(GeneralMatrix&);

void SortDescending(GeneralMatrix&);


void FFT(const ColumnVector&, const ColumnVector&,
   ColumnVector&, ColumnVector&);

void FFTI(const ColumnVector&, const ColumnVector&,
   ColumnVector&, ColumnVector&);

void RealFFT(const ColumnVector&, ColumnVector&, ColumnVector&);

void RealFFTI(const ColumnVector&, const ColumnVector&, ColumnVector&);


/********************** linear equation solving ****************************/

class LinearEquationSolver : public BaseMatrix
{
   GeneralMatrix* gm;
   int search(const BaseMatrix*) const;
   friend class BaseMatrix;
public:
   LinearEquationSolver(const BaseMatrix& bm);
   ~LinearEquationSolver();
   void CleanUp();
   GeneralMatrix* Evaluate(MatrixType);
   // probably should have an error message if MatrixType != UnSp
   NEW_DELETE(LinearEquationSolver)
};





#endif
