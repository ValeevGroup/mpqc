static char rcsid[] = "$Header$";

//------------------------------------------------------------------------
// Copyright (C) 1993: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//------------------------------------------------------------------------

#define WANT_STREAM
#define WANT_MATH
#include <math/newmat7/newmat.h>

static const char* separator = " ";

class PrintCounter
{
   int count;
   char* s;
public:
   ~PrintCounter();
   PrintCounter(char * sx) : count(0), s(sx) {}
   void operator++() { count++; }
};

   PrintCounter PCZ("Number of non-zero matrices = ");
   PrintCounter PCN("Number of matrices tested   = ");

PrintCounter::~PrintCounter()
{ cout << s << count << "\n"; }


void Print(const Matrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << (char*)X.Type() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<=nc; j++)  cout << X(i,j) << separator;
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const UpperTriangularMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << (char*)X.Type() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<i; j++) cout << separator;
      for (j=i; j<=nc; j++)  cout << X(i,j) << separator;
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const DiagonalMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << (char*)X.Type() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<i; j++) cout << separator;
      if (i<=nc) cout << X(i,i) << separator;
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const SymmetricMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << (char*)X.Type() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows(); int nc=X.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<i; j++) cout << X(j,i) << separator;
      for (j=i; j<=nc; j++)  cout << X(i,j) << separator;
      cout << "\n";
   }
   cout << flush; ++PCZ;
}

void Print(const LowerTriangularMatrix& X)
{
   ++PCN;
   cout << "\nMatrix type: " << (char*)X.Type() << " (";
   cout << X.Nrows() << ", ";
   cout << X.Ncols() << ")\n\n";
   if (X.IsZero()) { cout << "All elements are zero\n" << flush; return; }
   int nr=X.Nrows();
   for (int i=1; i<=nr; i++)
   {
      for (int j=1; j<=i; j++) cout << X(i,j) << separator;
      cout << "\n";
   }
   cout << flush; ++PCZ;
}


void Clean(Matrix& A, Real c)
{
   int nr = A.Nrows(); int nc = A.Ncols();
   for (int i=1; i<=nr; i++)
   {
      for ( int j=1; j<=nc; j++)
      { Real a = A(i,j); if ((a < c) && (a > -c)) A(i,j) = 0.0; }
   }
}

void Clean(DiagonalMatrix& A, Real c)
{
   int nr = A.Nrows();
   for (int i=1; i<=nr; i++)
   { Real a = A(i,i); if ((a < c) && (a > -c)) A(i,i) = 0.0; }
}


