//$$ sort.cxx                            Sorting

// Copyright (C) 1991,2,3: R B Davies

#define WANT_MATH

#include "include.h"

#include "newmatap.h"


/******************************** Shell sort ********************************/

void SortAscending(GeneralMatrix& GM)
{
   // from numerical recipies in C - Shell sort
   Tracer et("Sort-ascending");

   const double aln2i = 1.442695022; const double tiny = 1.0e-5;
   Real* gm = GM.Store(); int n = GM.Storage(); int m = n;
   int lognb2 = (int)(aln2i * log((double)n) + tiny);
   while (lognb2--)
   {
      m >>= 1;
      for (int j = m; j<n; j++)
      {
         Real* gmj = gm+j; int i = j-m; Real* gmi = gmj-m; Real t = *gmj;
         while (i>=0 && *gmi>t)  { *gmj = *gmi; gmj = gmi; gmi -= m; i -= m; }
         *gmj = t;
      }
   }
}

void SortDescending(GeneralMatrix& GM)
{
   // from numerical recipies in C - Shell sort
   Tracer et("Sort-descending");

   const double aln2i = 1.442695022; const double tiny = 1.0e-5;
   Real* gm = GM.Store(); int n = GM.Storage(); int m = n;
   int lognb2 = (int)(aln2i * log((double)n) + tiny);
   while (lognb2--)
   {
      m >>= 1;
      for (int j = m; j<n; j++)
      {
         Real* gmj = gm+j; int i = j-m; Real* gmi = gmj-m; Real t = *gmj;
         while (i>=0 && *gmi<t)  { *gmj = *gmi; gmj = gmi; gmi -= m; i -= m; }
         *gmj = t;
      }
   }
}

