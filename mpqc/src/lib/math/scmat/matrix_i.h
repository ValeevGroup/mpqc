
#ifndef _math_scmat_matrix_i_h
#define _math_scmat_matrix_i_h
#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/matrix.h>

// These are the inline candidates for the members defined in matrix.h.

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

///////////////////////////////////////////////////////////////////////////
// SCMatrixdouble inline candidates

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
// SymmSCMatrixdouble inline candidates

INLINE
SymmSCMatrixdouble::SymmSCMatrixdouble(SymmSCMatrix*a,int b,int c):
  matrix(a),i(b),j(c)
{
}
INLINE
SymmSCMatrixdouble::~SymmSCMatrixdouble()
{
}
INLINE double
SymmSCMatrixdouble::operator=(double a)
{
  matrix.set_element(i,j,a);
  return a;
}
INLINE
SymmSCMatrixdouble::operator double()
{
  return matrix.get_element(i,j);
}
INLINE double
SymmSCMatrixdouble::val()
{
  return matrix.get_element(i,j);
}

///////////////////////////////////////////////////////////////////////////
// DiagSCMatrixdouble inline candidates

INLINE
DiagSCMatrixdouble::DiagSCMatrixdouble(DiagSCMatrix*a,int b,int c):
  matrix(a),i(b),j(c)
{
}
INLINE
DiagSCMatrixdouble::~DiagSCMatrixdouble()
{
}
INLINE double
DiagSCMatrixdouble::operator=(double a)
{
  matrix.set_element(i,a);
  return a;
}
INLINE
DiagSCMatrixdouble::operator double()
{
  return matrix.get_element(i);
}
INLINE double
DiagSCMatrixdouble::val()
{
  return matrix.get_element(i);
}

///////////////////////////////////////////////////////////////////////////
// SCVectordouble inline candidates

INLINE
SCVectordouble::SCVectordouble(SCVector*a,int b):
  vector(a),i(b)
{
}
INLINE
SCVectordouble::~SCVectordouble()
{
}
INLINE double
SCVectordouble::operator=(double a)
{
  vector.set_element(i,a);
  return a;
}
INLINE
SCVectordouble::operator double()
{
  return vector.get_element(i);
}
INLINE double
SCVectordouble::val()
{
  return vector.get_element(i);
}

#undef INLINE

#endif
