#define WANT_STREAM

#include "include.h"

#include "newmat.h"

/**************************** test program ******************************/

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
      for (int j=1; j<=nc; j++)  cout << X(i,j) << "\t";
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
      for (int j=1; j<i; j++) cout << "\t";
      for (j=i; j<=nc; j++)  cout << X(i,j) << "\t";
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
      for (int j=1; j<i; j++) cout << "\t";
      if (i<=nc) cout << X(i,i) << "\t";
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
      for (int j=1; j<i; j++) cout << X(j,i) << "\t";
      for (j=i; j<=nc; j++)  cout << X(i,j) << "\t";
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
      for (int j=1; j<=i; j++) cout << X(i,j) << "\t";
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


void trymat1(); void trymat2(); void trymat3();
void trymat4(); void trymat5(); void trymat6();
void trymat7(); void trymat8(); void trymat9();
void trymata(); void trymatb(); void trymatc();
void trymatd(); void trymate(); void trymatf();
void trymatg(); void trymath(); void trymati();

main()
{
   Real* s1; Real* s2; Real* s3; Real* s4;
   cout << "\nBegin test\n";   // Forces cout to allocate memory at beginning
   { Matrix A1(40,200); s1 = A1.Store(); }
   { Matrix A1(1,1); s3 = A1.Store(); }
   {
      Tracer et("Matrix test program");

      Matrix A(25,150);
      {
	 int i;
	 RowVector A(8);
	 for (i=1;i<=7;i++) A(i)=0.0; A(8)=1.0;
	 Print(A);
      }
      cout << "\n";

      TestTypeAdd(); TestTypeMult(); TestTypeOrder();

      trymat1();
      trymat2();
      trymat3();
      trymat4();
      trymat5();
      trymat6();
      trymat7();
      trymat8();
      trymat9();
      trymata();
      trymatb();
      trymatc();
      trymatd();
// M.E. Colvin Made this change:
//      trymate();    Kills gcc 2.3.1 with signal 11A
      trymatf();
      trymatg();
      trymath();
      trymati();
   }
   { Matrix A1(40,200); s2 = A1.Store(); }
   cout << "\nChecking for lost memory: "
      << (unsigned long)s1 << " " << (unsigned long)s2 << " ";
   if (s1 != s2) cout << " - error\n"; else cout << " - ok\n\n";
   { Matrix A1(1,1); s4 = A1.Store(); }
   cout << "\nChecking for lost memory: "
      << (unsigned long)s3 << " " << (unsigned long)s4 << " ";
   if (s3 != s4) cout << " - error\n"; else cout << " - ok\n\n";
#ifdef DO_FREE_CHECK
   FreeCheck::Status();
#endif 
   return 0;
}
