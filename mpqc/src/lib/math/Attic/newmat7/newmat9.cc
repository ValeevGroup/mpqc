//$$ newmat9.cxx         Input and output

// Copyright (C) 1991,2,3: R B Davies


#define WANT_STREAM

#include "include.h"

#include "newmat.h"
#include "newmatrc.h"
#include "newmatio.h"

//#define REPORT { static ExeCounter ExeCount(__LINE__,9); ExeCount++; }

#define REPORT {}

ostream& operator<<(ostream& s, const BaseMatrix& X)
{
   GeneralMatrix* gm = ((BaseMatrix&)X).Evaluate(); operator<<(s, *gm);
   gm->tDelete(); return s;
}


ostream& operator<<(ostream& s, const GeneralMatrix& X)
{
   MatrixRow mr((GeneralMatrix*)&X, LoadOnEntry);
   int w = s.width();  int nr = X.Nrows();  long f = s.flags();
   s.setf(ios::fixed, ios::floatfield);
   for (int i=1; i<=nr; i++)
   {
      int skip = mr.skip;  int storage = mr.storage;
      Real* store = mr.store+skip;  skip *= w+1;
      while (skip--) s << " ";
      while (storage--) s << setw(w) << *store++ << " ";
      mr.Next();  s << "\n";
   }
   s << flush;  s.flags(f); return s;
}

