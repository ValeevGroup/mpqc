
#ifndef _math_scmat_matrix_i_h
#define _math_scmat_matrix_i_h

#include <math/scmat/abstract.h>

// These are the inline candidates for the members defined in matrix.h.

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

///////////////////////////////////////////////////////////////////////////
// SCdouble inline candidates

INLINE
SCMatrixdouble::SCMatrixdouble(SCMatrix*a,int b,int c):
  matrix(a),i(b),j(c)
{
}
INLINE
SCMatrixdouble::~SCMatrixdouble()
{
}
INLINE double
SCMatrixdouble::operator=(double a)
{
  matrix.set_element(i,j,a);
  return a;
}
INLINE
SCMatrixdouble::operator double()
{
  return matrix.get_element(i,j);
}
INLINE double
SCMatrixdouble::val()
{
  return matrix.get_element(i,j);
}

///////////////////////////////////////////////////////////////////////////
// RefSCMatrix inline candidates
// 
// INLINE RefSCMatrix :: RefSCMatrix() {}
// INLINE RefSCMatrix :: RefSCMatrix (RefSCMatrix & o):
//   RefSSSCMatrix (o)
// {}
// INLINE RefSCMatrix :: RefSCMatrix (SCMatrix * o):
//   RefSSSCMatrix (o)
// {}
// INLINE RefSCMatrix :: RefSCMatrix (RefDescribedClassBase&o):
//   RefSSSCMatrix (o)
// {}
// INLINE RefSCMatrix :: ~RefSCMatrix () {}
// INLINE RefSCMatrix& RefSCMatrix :: operator=(SCMatrix* cr)
// {
//   RefDCSCMatrix::operator=(cr);
//   return *this;
// }
// INLINE RefSCMatrix& RefSCMatrix :: operator=( RefDescribedClassBase & c)
// {
//   RefDCSCMatrix::operator=(c);
//   return *this;
// }
// INLINE RefSCMatrix& RefSCMatrix :: operator=( RefSCMatrix & c)
// {
//   RefDCSCMatrix::operator=(c);
//   return *this;
// }


#undef INLINE

#endif
