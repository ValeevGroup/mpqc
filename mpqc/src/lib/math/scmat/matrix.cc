
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/abstract.h>

/////////////////////////////////////////////////////////////////////////////
// SCElementOp member functions

SavableState_REF_def(SCElementOp);

#define CLASSNAME SCElementOp
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCElementOp::SCElementOp()
{
}

void *
SCElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCElementOp::~SCElementOp()
{
}

int
SCElementOp::has_collect()
{
  return 0;
}

void
SCElementOp::collect(RefSCElementOp&)
{
}

// If specializations of SCElementOp do not handle a particle
// block type, then these functions will be called and will
// set up an appropiate block iterator which specializations
// of SCElementOp must handle since it is pure virtual.

void
SCElementOp::process(SCMatrixRectBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixRectBlockIter(a);
  SCMatrixBlockIter&r=*i;
  process(r);
  // this causes a SCMatrixRectBlock::operator int() to be
  // called with this = 0x0 using gcc 2.5.6
  // process(*i,b);
  delete i;
}
void
SCElementOp::process(SCMatrixLTriBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixLTriBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process(SCMatrixDiagBlock* a)
{
  SCMatrixBlockIter*i = new SCMatrixDiagBlockIter(a);
  process(*i);
  delete i;
}
void
SCElementOp::process(SCVectorSimpleBlock* a)
{
  SCMatrixBlockIter*i = new SCVectorSimpleBlockIter(a);
  process(*i);
  delete i;
}

/////////////////////////////////////////////////////////////////////////////
// SCRectElementOp member functions

#define CLASSNAME SCRectElementOp
#define PARENTS virtual public SCElementOp
#include <util/state/statei.h>
#include <util/class/classia.h>

SavableState_REF_def(SCRectElementOp);

void *
SCRectElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCRectElementOp::SCRectElementOp()
{
}

SCRectElementOp::~SCRectElementOp()
{
}

/////////////////////////////////////////////////////////////////////////////
// SCSymmElementOp member functions

#define CLASSNAME SCSymmElementOp
#define PARENTS virtual public SCElementOp
#include <util/state/statei.h>
#include <util/class/classia.h>

SavableState_REF_def(SCSymmElementOp);

void *
SCSymmElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCSymmElementOp::SCSymmElementOp()
{
}

SCSymmElementOp::~SCSymmElementOp()
{
}

/////////////////////////////////////////////////////////////////////////////
// SCDiagElementOp member functions

#define CLASSNAME SCDiagElementOp
#define PARENTS virtual public SCElementOp
#include <util/state/statei.h>
#include <util/class/classia.h>

SavableState_REF_def(SCDiagElementOp);

void *
SCDiagElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCDiagElementOp::SCDiagElementOp()
{
}

SCDiagElementOp::~SCDiagElementOp()
{
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorElementOp member functions

#define CLASSNAME SCVectorElementOp
#define PARENTS virtual public SCElementOp
#include <util/state/statei.h>
#include <util/class/classia.h>

SavableState_REF_def(SCVectorElementOp);

void *
SCVectorElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCVectorElementOp::SCVectorElementOp()
{
}

SCVectorElementOp::~SCVectorElementOp()
{
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrix reference base class member functions

SavableState_named_REF_def(RefSSSCDimension,SCDimension);
SavableState_named_REF_def(RefSSSCVector,SCVector);
SavableState_named_REF_def(RefSSSCMatrix,SCMatrix);
SavableState_named_REF_def(RefSSSymmSCMatrix,SymmSCMatrix);
SavableState_named_REF_def(RefSSDiagSCMatrix,DiagSCMatrix);

/////////////////////////////////////////////////////////////////////////////
// SCDimension reference member functions

RefSCDimension::RefSCDimension() {}
             
RefSCDimension::RefSCDimension (RefSCDimension & o): RefSSSCDimension (o) {}
             
RefSCDimension::RefSCDimension (SCDimension * o): RefSSSCDimension (o) {}
             
RefSCDimension::RefSCDimension (RefDescribedClassBase&o): RefSSSCDimension (o) {}

RefSCDimension::~RefSCDimension () {}

RefSCDimension&
RefSCDimension::operator=(SCDimension* cr)
{
  RefSSSCDimension::operator=(cr);
  return *this;
}

RefSCDimension&
RefSCDimension::operator=( RefDescribedClassBase & c)
{
  RefSSSCDimension::operator=(c);
  return *this;
}

RefSCDimension&
RefSCDimension::operator=( RefSCDimension & c)
{
  RefSSSCDimension::operator=(c);
  return *this;
}

int
RefSCDimension::n()
{
  if (null()) return 0;
  return pointer()->n();
}

RefSCDimension::operator int()
{
  if (null()) return 0;
  return pointer()->n();
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrix reference member functions
RefSCMatrix::RefSCMatrix() {}
             
RefSCMatrix::RefSCMatrix (RefSCMatrix & o): RefSSSCMatrix (o) {}
             
RefSCMatrix::RefSCMatrix (StateIn & o): RefSSSCMatrix (o) {}
             
RefSCMatrix::RefSCMatrix (SCMatrix * o): RefSSSCMatrix (o) {}
             
// RefSCMatrix::RefSCMatrix (RefDescribedClassBase&o): RefSSSCMatrix (o) {}

RefSCMatrix::~RefSCMatrix () {}

RefSCMatrix&
RefSCMatrix::operator=(SCMatrix* cr)
{
  RefSSSCMatrix::operator=(cr);
  return *this;
}

// RefSCMatrix&
// RefSCMatrix::operator=( RefDescribedClassBase & c)
// {
//   RefSSSCMatrix::operator=(c);
//   return *this;
// }

RefSCMatrix&
RefSCMatrix::operator=( RefSCMatrix & c)
{
  RefSSSCMatrix::operator=(c);
  return *this;
}

RefSCMatrix::RefSCMatrix(RefSCDimension&a,RefSCDimension&b)
{
  a.require_nonnull();
  b.require_nonnull();
  assign_pointer(a->create_matrix(b.pointer()));
}

void
RefSCMatrix::set_element(int i, int j, double a)
{
  require_nonnull();
  pointer()->set_element(i,j,a);
}

double
RefSCMatrix::get_element(int i, int j)
{
  require_nonnull();
  return pointer()->get_element(i,j);
}

RefSCVector
RefSCMatrix::operator*(RefSCVector&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector r = rowdim()->create_vector();
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator*(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = rowdim()->create_matrix(a->coldim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator*(RefSymmSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = rowdim()->create_matrix(a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator*(RefDiagSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = rowdim()->create_matrix(a->dim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCMatrix
RefSCMatrix::operator+(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix ret(rowdim(),coldim());
  
  ret->assign(pointer());
  ret->accumulate(a);

  return ret;
}

RefSCMatrix
RefSCMatrix::operator-(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix ret(rowdim(),coldim());
  
  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

RefSCMatrix
RefSCMatrix::t()
{
  require_nonnull();
  
  RefSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->transpose_this();
  return ret;
}

RefSCMatrix
RefSCMatrix::i()
{
  require_nonnull();
  
  RefSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->invert_this();
  return ret;
}

int
RefSCMatrix::nrow()
{
  if (null()) return 0;
  else return pointer()->nrow();
}

int
RefSCMatrix::ncol()
{
  if (null()) return 0;
  else return pointer()->ncol();
}

RefSCDimension
RefSCMatrix::rowdim()
{
  if (null()) return 0;
  else return pointer()->rowdim();
}

RefSCDimension
RefSCMatrix::coldim()
{
  if (null()) return 0;
  else return pointer()->coldim();
}

SCMatrixdouble
RefSCMatrix::operator()(int i,int j)
{
  return SCMatrixdouble(pointer(),i,j);
}

RefSCMatrix
RefSCMatrix::clone()
{
  RefSCMatrix r = rowdim()->create_matrix(coldim());
  return r;
}

void
RefSCMatrix::accumulate_product(RefSCMatrix&a,RefSCMatrix&b)
{
  require_nonnull();
  pointer()->accumulate_product(a.pointer(),b.pointer());
}

RefSCMatrix
RefSCMatrix::copy()
{
  if (null()) return 0;
  RefSCMatrix v = rowdim()->create_matrix(coldim());
  v.assign(*this);
  return v;
}

void
RefSCMatrix::assign(RefSCMatrix&a)
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefSCMatrix::assign(const double*v)
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSCMatrix::assign(const double**v)
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSCMatrix::convert(double*v)
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSCMatrix::convert(double**v)
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSCMatrix::scale(double a)
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefSCMatrix::assign(double a)
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefSCMatrix::accumulate(RefSCMatrix&a)
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSCMatrix::element_op(RefSCRectElementOp&op)
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefSCMatrix::print(ostream& out)
{
  print(0,out);
}

void
RefSCMatrix::print(const char*title,ostream&out, int precision)
{
  require_nonnull();
  pointer()->print(title,out,precision);
}

RefSCMatrix
RefSCMatrix::operator *(double a)
{
  RefSCMatrix r(copy());
  r.scale(a);
  return r;
}

RefSCMatrix
operator *(double a, RefSCMatrix& v)
{
  return v*a;
}

void
RefSCMatrix::accumulate_outer_product(RefSCVector& v1,RefSCVector&v2)
{
  require_nonnull();
  pointer()->accumulate_outer_product(v1.pointer(),v2.pointer());
}

///////////////////////////////////////////////////////////////////
// RefSymmSCMatrix members

RefSymmSCMatrix::RefSymmSCMatrix()
{
}
             
RefSymmSCMatrix::RefSymmSCMatrix (StateIn & o):
  RefSSSymmSCMatrix (o)
{
}
             
RefSymmSCMatrix::RefSymmSCMatrix (RefSymmSCMatrix & o):
  RefSSSymmSCMatrix (o)
{
}
             
RefSymmSCMatrix::RefSymmSCMatrix (SymmSCMatrix * o):
  RefSSSymmSCMatrix (o)
{
}
             
// RefSymmSCMatrix::RefSymmSCMatrix (RefDescribedClassBase&o):
//   RefSSSymmSCMatrix (o)
// {
// }

RefSymmSCMatrix::~RefSymmSCMatrix ()
{
}

RefSymmSCMatrix&
RefSymmSCMatrix::operator=(SymmSCMatrix* cr)
{
  RefSSSymmSCMatrix::operator=(cr);
  return *this;
}

// RefSymmSCMatrix&
// RefSymmSCMatrix::operator=( RefDescribedClassBase & c)
// {
//   RefSSSymmSCMatrix::operator=(c);
//   return *this;
// }

RefSymmSCMatrix&
RefSymmSCMatrix::operator=( RefSymmSCMatrix & c)
{
  RefSSSymmSCMatrix::operator=(c);
  return *this;
}

RefSymmSCMatrix::RefSymmSCMatrix(RefSCDimension&a)
{
  a.require_nonnull();
  assign_pointer(a->create_symmmatrix());
}

void
RefSymmSCMatrix::set_element(int i, int j, double a)
{
  require_nonnull();
  pointer()->set_element(i,j,a);
}

double
RefSymmSCMatrix::get_element(int i, int j)
{
  require_nonnull();
  return pointer()->get_element(i,j);
}

void
RefSymmSCMatrix::accumulate_symmetric_product(RefSCMatrix& a)
{
  require_nonnull();
  pointer()->accumulate_symmetric_product(a.pointer());
}

void
RefSymmSCMatrix::accumulate_symmetric_sum(RefSCMatrix& a)
{
  require_nonnull();
  pointer()->accumulate_symmetric_sum(a.pointer());
}

void
RefSymmSCMatrix::accumulate_transform(RefSCMatrix& a,RefSymmSCMatrix&b)
{
  require_nonnull();
  pointer()->accumulate_transform(a.pointer(),b.pointer());
}

RefSymmSCMatrix
RefSymmSCMatrix::operator+(RefSymmSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSymmSCMatrix ret(dim());
  
  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::operator-(RefSymmSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSymmSCMatrix ret(dim());
  
  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::i()
{
  require_nonnull();
  
  RefSymmSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->invert_this();
  return ret;
}

int
RefSymmSCMatrix::n()
{
  if (null()) return 0;
  else return pointer()->dim()->n();
}

RefSCDimension
RefSymmSCMatrix::dim()
{
  if (null()) return 0;
  else return pointer()->dim();
}

SymmSCMatrixdouble
RefSymmSCMatrix::operator()(int i,int j)
{
  return SymmSCMatrixdouble(pointer(),i,j);
}

RefSymmSCMatrix
RefSymmSCMatrix::clone()
{
  RefSymmSCMatrix r = dim()->create_symmmatrix();
  return r;
}

RefDiagSCMatrix
RefSymmSCMatrix::eigvals()
{
  if (null()) return 0;
  RefDiagSCMatrix vals = dim()->create_diagmatrix();
  RefSCMatrix vecs = dim()->create_matrix(dim());
  diagonalize(vals,vecs);
  return vals;
}

RefSCMatrix
RefSymmSCMatrix::eigvecs()
{
  if (null()) return 0;
  RefDiagSCMatrix vals = dim()->create_diagmatrix();
  RefSCMatrix vecs = dim()->create_matrix(dim());
  diagonalize(vals,vecs);
  return vecs;
}

void
RefSymmSCMatrix::diagonalize(RefDiagSCMatrix& vals, RefSCMatrix& vecs)
{
  require_nonnull();
  pointer()->diagonalize(vals.pointer(),vecs.pointer());
}

RefSymmSCMatrix
RefSymmSCMatrix::copy()
{
  if (null()) return 0;
  RefSymmSCMatrix v = dim()->create_symmmatrix();
  v.assign(*this);
  return v;
}

void
RefSymmSCMatrix::assign(RefSymmSCMatrix&a)
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefSymmSCMatrix::assign(const double*v)
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSymmSCMatrix::assign(const double**v)
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSymmSCMatrix::convert(double*v)
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSymmSCMatrix::convert(double**v)
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSymmSCMatrix::scale(double a)
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefSymmSCMatrix::assign(double a)
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefSymmSCMatrix::accumulate(RefSymmSCMatrix&a)
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSymmSCMatrix::element_op(RefSCSymmElementOp&op)
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefSymmSCMatrix::print(ostream& out)
{
  print(0,out);
}

void
RefSymmSCMatrix::print(const char*title,ostream&out, int precision)
{
  require_nonnull();
  pointer()->print(title,out,precision);
}

RefSCMatrix
RefSymmSCMatrix::operator*(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = dim()->create_matrix(a->coldim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSCVector
RefSymmSCMatrix::operator*(RefSCVector&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector r = dim()->create_vector();
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefSymmSCMatrix
RefSymmSCMatrix::operator *(double a)
{
  RefSymmSCMatrix r(copy());
  r.scale(a);
  return r;
}

RefSymmSCMatrix
operator *(double a, RefSymmSCMatrix& v)
{
  return v*a;
}

void
RefSymmSCMatrix::accumulate_symmetric_outer_product(RefSCVector&v)
{
  require_nonnull();
  pointer()->accumulate_symmetric_outer_product(v.pointer());
}

double
RefSymmSCMatrix::scalar_product(RefSCVector&v)
{
  if (null()) return 0.0;
  return pointer()->scalar_product(v.pointer());
}

///////////////////////////////////////////////////////////////////
// RefDiagSCMatrix members

RefDiagSCMatrix::RefDiagSCMatrix()
{
}
             
RefDiagSCMatrix::RefDiagSCMatrix (StateIn & o):
  RefSSDiagSCMatrix (o)
{
}
             
RefDiagSCMatrix::RefDiagSCMatrix (RefDiagSCMatrix & o):
  RefSSDiagSCMatrix (o)
{
}
             
RefDiagSCMatrix::RefDiagSCMatrix (DiagSCMatrix * o):
  RefSSDiagSCMatrix (o)
{
}
             
// RefDiagSCMatrix::RefDiagSCMatrix (RefDescribedClassBase&o):
//   RefSSDiagSCMatrix (o)
// {
// }

RefDiagSCMatrix::~RefDiagSCMatrix ()
{
}

RefDiagSCMatrix&
RefDiagSCMatrix::operator=(DiagSCMatrix* cr)
{
  RefSSDiagSCMatrix::operator=(cr);
  return *this;
}

// RefDiagSCMatrix&
// RefDiagSCMatrix::operator=( RefDescribedClassBase & c)
// {
//   RefSSDiagSCMatrix::operator=(c);
//   return *this;
// }

RefDiagSCMatrix&
RefDiagSCMatrix::operator=( RefDiagSCMatrix & c)
{
  RefSSDiagSCMatrix::operator=(c);
  return *this;
}

RefDiagSCMatrix::RefDiagSCMatrix(RefSCDimension&a)
{
  a.require_nonnull();
  assign_pointer(a->create_diagmatrix());
}

void
RefDiagSCMatrix::set_element(int i, double a)
{
  require_nonnull();
  pointer()->set_element(i,a);
}

double
RefDiagSCMatrix::get_element(int i)
{
  require_nonnull();
  return pointer()->get_element(i);
}

RefSCMatrix
RefDiagSCMatrix::operator*(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix r = dim()->create_matrix(a->coldim());
  r->assign(0.0);
  r->accumulate_product(pointer(),a.pointer());
  return r;
}

RefDiagSCMatrix
RefDiagSCMatrix::operator+(RefDiagSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix ret(dim());
  
  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefDiagSCMatrix
RefDiagSCMatrix::operator-(RefDiagSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix ret(dim());
  
  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

RefDiagSCMatrix
RefDiagSCMatrix::i()
{
  require_nonnull();
  
  RefDiagSCMatrix ret;
  ret = clone();
  ret->assign(pointer());
  ret->invert_this();
  return ret;
}

int
RefDiagSCMatrix::n()
{
  if (null()) return 0;
  else return pointer()->dim()->n();
}

RefSCDimension
RefDiagSCMatrix::dim()
{
  if (null()) return 0;
  else return pointer()->dim();
}

DiagSCMatrixdouble
RefDiagSCMatrix::operator()(int i)
{
  return DiagSCMatrixdouble(pointer(),i,i);
}

RefDiagSCMatrix
RefDiagSCMatrix::clone()
{
  RefDiagSCMatrix r = dim()->create_diagmatrix();
  return r;
}

RefDiagSCMatrix
RefDiagSCMatrix::copy()
{
  if (null()) return 0;
  RefDiagSCMatrix v = dim()->create_diagmatrix();
  v.assign(*this);
  return v;
}

void
RefDiagSCMatrix::assign(RefDiagSCMatrix&a)
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefDiagSCMatrix::assign(const double*v)
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefDiagSCMatrix::convert(double*v)
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefDiagSCMatrix::scale(double a)
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefDiagSCMatrix::assign(double a)
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefDiagSCMatrix::accumulate(RefDiagSCMatrix&a)
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefDiagSCMatrix::element_op(RefSCDiagElementOp&op)
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefDiagSCMatrix::print(ostream& out)
{
  print(0,out);
}

void
RefDiagSCMatrix::print(const char*title,ostream&out, int precision)
{
  require_nonnull();
  pointer()->print(title,out,precision);
}

RefDiagSCMatrix
RefDiagSCMatrix::operator *(double a)
{
  RefDiagSCMatrix r(copy());
  r.scale(a);
  return r;
}

RefDiagSCMatrix
operator *(double a, RefDiagSCMatrix& v)
{
  return v*a;
}

///////////////////////////////////////////////////////////////////
// RefSCVector members

ARRAY_def(RefSCVector);
SET_def(RefSCVector);

RefSCVector::RefSCVector()
{
}
             
RefSCVector::RefSCVector (StateIn & o):
  RefSSSCVector (o)
{
}
             
RefSCVector::RefSCVector (RefSCVector & o):
  RefSSSCVector (o)
{
}
             
RefSCVector::RefSCVector (SCVector * o):
  RefSSSCVector (o)
{
}

// RefSCVector::RefSCVector (RefDescribedClassBase&o):
//   RefSSSCVector (o)
// {
// }

RefSCVector::~RefSCVector ()
{
}

RefSCVector&
RefSCVector::operator=(SCVector* cr)
{
  RefSSSCVector::operator=(cr);
  return *this;
}

// RefSCVector&
// RefSCVector::operator=( RefDescribedClassBase & c)
// {
//   RefSSSCVector::operator=(c);
//   return *this;
// }

RefSCVector&
RefSCVector::operator=( RefSCVector & c)
{
  RefSSSCVector::operator=(c);
  return *this;
}

RefSCVector::RefSCVector(RefSCDimension&a)
{
  a.require_nonnull();
  assign_pointer(a->create_vector());
}

void
RefSCVector::set_element(int i, double a)
{
  require_nonnull();
  pointer()->set_element(i,a);
}

double
RefSCVector::get_element(int i)
{
  require_nonnull();
  return pointer()->get_element(i);
}

RefSCVector
RefSCVector::operator+(RefSCVector&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector ret(dim());
  
  ret->assign(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefSCVector
RefSCVector::operator-(RefSCVector&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCVector ret(dim());
  
  ret->assign(a.pointer());
  ret->scale(-1.0);
  ret->accumulate(pointer());

  return ret;
}

int
RefSCVector::n()
{
  if (null()) return 0;
  else return pointer()->dim()->n();
}

RefSCDimension
RefSCVector::dim()
{
  if (null()) return 0;
  else return pointer()->dim();
}

SCVectordouble
RefSCVector::operator()(int i)
{
  return SCVectordouble(pointer(),i);
}

SCVectordouble
RefSCVector::operator[](int i)
{
  return SCVectordouble(pointer(),i);
}

RefSCVector
RefSCVector::clone()
{
  RefSCVector r = dim()->create_vector();
  return r;
}

RefSCVector
RefSCVector::copy()
{
  if (null()) return 0;
  RefSCVector v = dim()->create_vector();
  v.assign(*this);
  return v;
}

double
RefSCVector::dot(RefSCVector&a)
{
  require_nonnull();
  return pointer()->scalar_product(a.pointer());
}

double
RefSCVector::scalar_product(RefSCVector&a)
{
  require_nonnull();
  return pointer()->scalar_product(a.pointer());
}

void
RefSCVector::assign(RefSCVector&a)
{
  require_nonnull();
  pointer()->assign(a.pointer());
}

void
RefSCVector::assign(const double*v)
{
  require_nonnull();
  pointer()->assign(v);
}

void
RefSCVector::convert(double*v)
{
  require_nonnull();
  pointer()->convert(v);
}

void
RefSCVector::scale(double a)
{
  require_nonnull();
  pointer()->scale(a);
}

void
RefSCVector::assign(double a)
{
  require_nonnull();
  pointer()->assign(a);
}

void
RefSCVector::accumulate(RefSCVector&a)
{
  require_nonnull();
  pointer()->accumulate(a.pointer());
}

void
RefSCVector::element_op(RefSCVectorElementOp&op)
{
  if (nonnull()) pointer()->element_op(op);
}

void
RefSCVector::print(ostream& out)
{
  print(0,out);
}

void
RefSCVector::print(const char*title,ostream&out, int precision)
{
  require_nonnull();
  pointer()->print(title,out,precision);
}

RefSCVector
RefSCVector::operator *(double a)
{
  RefSCVector r(copy());
  r.scale(a);
  return r;
}

RefSCVector
operator *(double a, RefSCVector& v)
{
  return v*a;
}

void
RefSCVector::normalize()
{
  require_nonnull();
  pointer()->normalize();
}

RefSymmSCMatrix
RefSCVector::symmetric_outer_product()
{
  RefSymmSCMatrix result(dim());
  result.assign(0.0);
  result.accumulate_symmetric_outer_product(pointer());
  return result;
}

RefSCMatrix
RefSCVector::outer_product(RefSCVector&v)
{
  RefSCMatrix result(dim(),v.dim());
  result.assign(0.0);
  result.accumulate_outer_product(*this,v);
  return result;
}

double
RefSCVector::maxabs()
{
  if (null()) return 0.0;
  return pointer()->maxabs();
}

