//$$ newmatex.cxx                    Exception handler

// Copyright (C) 1992: R B Davies

#define WANT_STREAM                  // include.h will get stream fns

#include "include.h"                 // include standard files
#include "newmat.h"


// action = -1    print message and exit(1)
//           0    no message if handler available
//           1    print message and use handler



int SpaceException::action = 1;
int DataException::action = 1;
int ConvergenceException::action = 1;
int ProgramException::action = 1;
int InternalException::action = 1;


static inline iabs(int i) { return i >= 0 ? i : -i; }


MatrixDetails::MatrixDetails(const GeneralMatrix& A)
   : type(A.Type()), nrows(A.Nrows()), ncols(A.Ncols())
{ MatrixBandWidth bw = A.BandWidth(); ubw = bw.upper; lbw = bw.lower; }

void MatrixDetails::PrintOut()
{
   cout << "MatrixType = " << (char*)type;
   cout << "  # Rows = " << nrows;
   cout << "; # Cols = " << ncols;
   if (lbw >=0) cout << "; lower BW = " << lbw;
   if (ubw >=0) cout << "; upper BW = " << ubw;
   cout << "\n";
}



SpaceException::SpaceException() : Exception(iabs(action))
{
   if (action) cout << "Out of space on heap\n";
   if (action < 0) exit(1);
}

MatrixException::MatrixException(int action) : Exception(iabs(action))
{ if (action) cout << "The exception is from newmat.\n"; }

MatrixException::MatrixException(int action, const GeneralMatrix& A)
   : Exception(iabs(action))
{
   if (action)
   {
      cout << "The exception is from newmat: details of matrix follow:\n";
      MatrixDetails(A).PrintOut();
   }
}

MatrixException::MatrixException(int action, const GeneralMatrix& A,
   const GeneralMatrix& B) : Exception(iabs(action))
{
   if (action)
   {
      cout << "The exception is from newmat: details of matrices follow:\n";
      MatrixDetails(A).PrintOut();
      MatrixDetails(B).PrintOut();
   }
}

DataException::DataException(const GeneralMatrix& A)
   : MatrixException(action, A) {}

NPDException::NPDException(const GeneralMatrix& A)
   : DataException(A)
{
   if (action) cout << "The matrix is not positive definite\n\n";
   if (action < 0) exit(1);
}

SingularException::SingularException(const GeneralMatrix& A)
   : DataException(A)
{
   if (action) cout << "The matrix is singular\n\n";
   if (action < 0) exit(1);
}

ConvergenceException::ConvergenceException(const GeneralMatrix& A)
   : MatrixException(action,A)
{
   if (action) cout << "Process fails to converge\n\n";
   if (action < 0) exit(1);
}

ProgramException::ProgramException(char* c) : MatrixException(action)
{
   if (action) cout << c << "\n\n";
   if (action < 0) exit(1);
}

ProgramException::ProgramException(char* c, const GeneralMatrix& A)
   : MatrixException(action,A)
{
   if (action) cout << c << "\n\n";
   if (action < 0) exit(1);
}

ProgramException::ProgramException(char* c, const GeneralMatrix& A,
   const GeneralMatrix& B) : MatrixException(action,A,B)
{
   if (action) cout << c << "\n\n";
   if (action < 0) exit(1);
}

ProgramException::ProgramException(const GeneralMatrix& A)
   : MatrixException(action, A) {}

ProgramException::ProgramException() : MatrixException(action) {}

VectorException::VectorException() : ProgramException()
{
   if (action) cout << "Cannot convert matrix to vector\n\n";
   if (action < 0) exit(1);
}

VectorException::VectorException(const GeneralMatrix& A)
   : ProgramException(A)
{
   if (action) cout << "Cannot convert matrix to vector\n\n";
   if (action < 0) exit(1);
}

NotSquareException::NotSquareException(const GeneralMatrix& A)
   : ProgramException(A)
{
   if (action) cout << "Matrix is not square\n\n";
   if (action < 0) exit(1);
}

SubMatrixDimensionException::SubMatrixDimensionException()
   : ProgramException()
{
   if (action) cout << "Incompatible submatrix dimension\n\n";
   if (action < 0) exit(1);
}

IncompatibleDimensionsException::IncompatibleDimensionsException()
   : ProgramException()
{
   if (action) cout << "Incompatible dimensions\n\n";
   if (action < 0) exit(1);
}

NotDefinedException::NotDefinedException(char* op, char* matrix)
   : ProgramException()
{
   if (action)
      cout << "Operation " << op << " not defined for " << matrix << "\n\n";
   if (action < 0) exit(1);
}

CannotBuildException::CannotBuildException(char* matrix)
   : ProgramException()
{
   if (action)
      cout << "Cannot build matrix type " << matrix << "\n\n";
   if (action < 0) exit(1);
}

IndexException::IndexException(int i, const GeneralMatrix& A)
   : ProgramException(A)
{
   if (action)
      { cout << "Index error: requested index = " << i << "\n\n"; }
   if (action < 0) exit(1);
}

IndexException::IndexException(int i, int j, const GeneralMatrix& A)
   : ProgramException(A)
{
   if (action)
   {
      cout << "Index error: requested indices = " << i << ", " << j << "\n\n";
   }
   if (action < 0) exit(1);
}


IndexException::IndexException(int i, const GeneralMatrix& A, Boolean)
   : ProgramException(A)
{
   if (action)
      { cout << "Element error: requested index (wrt 0) = " << i << "\n\n"; }
   if (action < 0) exit(1);
}

IndexException::IndexException(int i, int j, const GeneralMatrix& A, Boolean)
   : ProgramException(A)
{
   if (action)
   {
      cout << "Element error: requested indices (wrt 0) = " 
         << i << ", " << j << "\n\n";
   }
   if (action < 0) exit(1);
}

InternalException::InternalException(char* c) : MatrixException(action)
{
   if (action) cout << c << "\n\n";
   if (action < 0) exit(1);
}




/************************* ExeCounter functions *****************************/



int ExeCounter::nreports = 0;

ExeCounter::ExeCounter(int xl, int xf) : line(xl), fileid(xf), nexe(0) {}

ExeCounter::~ExeCounter()
{
   nreports++;
   cout << nreports << "  " << fileid << "  " << line << "  " << nexe << "\n";
}


 
/**************************** error handler *******************************/

void MatrixErrorNoSpace(void* v) { if (!v) Throw(SpaceException()); }
// throw exception if v is null



/************************* test type manipulation **************************/




// These functions may cause problems for Glockenspiel 2.0c; they are used
// only for testing so you can delete them


void TestTypeAdd()
{
   MatrixType list[] = { MatrixType::UT,
                         MatrixType::LT,
                         MatrixType::Rt,
                         MatrixType::Sm,
			 MatrixType::Dg,
                         MatrixType::BM,
                         MatrixType::UB,
			 MatrixType::LB,
			 MatrixType::SB };

   cout << "+     ";
   for (int i=0; i<MatrixType::nTypes(); i++) cout << (char*)list[i] << " ";
   cout << "\n";
   for (i=0; i<MatrixType::nTypes(); i++)
   {
      cout << (char*)(list[i]) << " ";
      for (int j=0; j<MatrixType::nTypes(); j++)
	 cout << (char*)(list[j]+list[i]) << " ";
      cout << "\n";
   }
   cout << "\n";
}

void TestTypeMult()
{
   MatrixType list[] = { MatrixType::UT,
                         MatrixType::LT,
                         MatrixType::Rt,
                         MatrixType::Sm,
                         MatrixType::Dg,
                         MatrixType::BM,
                         MatrixType::UB,
                         MatrixType::LB,
                         MatrixType::SB };
   cout << "*     ";
   for (int i=0; i<MatrixType::nTypes(); i++)
      cout << (char*)list[i] << " ";
   cout << "\n";
   for (i=0; i<MatrixType::nTypes(); i++)
   {
      cout << (char*)list[i] << " ";
      for (int j=0; j<MatrixType::nTypes(); j++)
	 cout << (char*)(list[j]*list[i]) << " ";
      cout << "\n";
   }
   cout << "\n";
}

void TestTypeOrder()
{
   MatrixType list[] = { MatrixType::UT,
                         MatrixType::LT,
                         MatrixType::Rt,
                         MatrixType::Sm,
                         MatrixType::Dg,
                         MatrixType::BM,
                         MatrixType::UB,
                         MatrixType::LB,
                         MatrixType::SB };
   cout << ">=    ";
   for (int i = 0; i<MatrixType::nTypes(); i++)
      cout << (char*)list[i] << " ";
   cout << "\n";
   for (i=0; i<MatrixType::nTypes(); i++)
   {
      cout << (char*)list[i] << " ";
      for (int j=0; j<MatrixType::nTypes(); j++)
	 cout << ((list[j]>=list[i]) ? "Yes   " : "No    ");
      cout << "\n";
   }
   cout << "\n";
}


/************************* miscellanous errors ***************************/


void CroutMatrix::GetRow(MatrixRowCol&)
   { Throw(NotDefinedException("GetRow","Crout")); }
void CroutMatrix::GetCol(MatrixRowCol&)
   { Throw(NotDefinedException("GetCol","Crout")); }
void CroutMatrix::operator=(const BaseMatrix&)
   { Throw(NotDefinedException("=","Crout")); }
void BandLUMatrix::GetRow(MatrixRowCol&)
   { Throw(NotDefinedException("GetRow","BandLUMatrix")); }
void BandLUMatrix::GetCol(MatrixRowCol&)
   { Throw(NotDefinedException("GetCol","BandLUMatrix")); }
void BandLUMatrix::operator=(const BaseMatrix&)
   { Throw(NotDefinedException("=","BandLUMatrix")); }
#ifdef TEMPS_DESTROYED_QUICKLY
   ReturnMatrixX::ReturnMatrixX(const ReturnMatrixX& tm)
     : gm(tm.gm)
   {
     Throw(ProgramException("ReturnMatrixX error"));
   }
#endif

