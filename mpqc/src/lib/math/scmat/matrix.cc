
#include <math/scmat/blkiter.h>
#include <math/scmat/abstract.h>

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlock member functions

#define CLASSNAME SCMatrixBlock
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCMatrixBlock::SCMatrixBlock()
{
}

void *
SCMatrixBlock::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCMatrixBlock::~SCMatrixBlock()
{
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixRectBlock member functions

SavableState_REF_def(SCMatrixRectBlock);

#define CLASSNAME SCMatrixRectBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCMatrixRectBlock::SCMatrixRectBlock(int is, int ie, int js, int je):
  istart(is),
  jstart(js),
  iend(ie),
  jend(je)
{
  data = new double[(ie-is)*(je-js)];
}

SCMatrixRectBlock::SCMatrixRectBlock(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(jend);
  s.get(data);
}

void
SCMatrixRectBlock::save_data_state(StateOut&s)
{
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(jend);
  s.put(data,(iend-istart)*(jend-jstart));
}

void *
SCMatrixRectBlock::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCMatrixBlock::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCMatrixRectBlock::~SCMatrixRectBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixLTriBlock member functions

SavableState_REF_def(SCMatrixLTriBlock);

#define CLASSNAME SCMatrixLTriBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCMatrixLTriBlock::SCMatrixLTriBlock(int s,int e):
  start(s),
  end(e)
{
  data = new double[((e-s)*(e-s+1))/2];
}

SCMatrixLTriBlock::SCMatrixLTriBlock(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(start);
  s.get(end);
  s.get(data);
}

void
SCMatrixLTriBlock::save_data_state(StateOut&s)
{
  s.put(start);
  s.put(end);
  s.put(data,((end-start)*(end-start+1))/2);
}

void *
SCMatrixLTriBlock::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCMatrixBlock::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCMatrixLTriBlock::~SCMatrixLTriBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixDiagBlock member functions

SavableState_REF_def(SCMatrixDiagBlock);

#define CLASSNAME SCMatrixDiagBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCMatrixDiagBlock::SCMatrixDiagBlock(int s, int e):
  istart(s),
  iend(e),
  jstart(s)
{
  data = new double[e-s];
}

SCMatrixDiagBlock::SCMatrixDiagBlock(int is, int ie,int js):
  istart(is),
  iend(ie),
  jstart(js)
{
  data = new double[ie-is];
}

SCMatrixDiagBlock::SCMatrixDiagBlock(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(istart);
  s.get(jstart);
  s.get(iend);
  s.get(data);
}

void
SCMatrixDiagBlock::save_data_state(StateOut&s)
{
  s.put(istart);
  s.put(jstart);
  s.put(iend);
  s.put(data,iend-istart);
}

void *
SCMatrixDiagBlock::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCMatrixBlock::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCMatrixDiagBlock::~SCMatrixDiagBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorSimpleBlock member functions

SavableState_REF_def(SCVectorSimpleBlock);

#define CLASSNAME SCVectorSimpleBlock
#define PARENTS public SCMatrixBlock
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCVectorSimpleBlock::SCVectorSimpleBlock(int s, int e):
  istart(s),
  iend(e)
{
  data = new double[e-s];
}

SCVectorSimpleBlock::SCVectorSimpleBlock(int is, int ie):
  istart(is),
  iend(ie)
{
  data = new double[ie-is];
}

SCVectorSimpleBlock::SCVectorSimpleBlock(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(istart);
  s.get(iend);
  s.get(data);
}

void
SCVectorSimpleBlock::save_data_state(StateOut&s)
{
  s.put(istart);
  s.put(iend);
  s.put(data,iend-istart);
}

void *
SCVectorSimpleBlock::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCMatrixBlock::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCVectorSimpleBlock::~SCVectorSimpleBlock()
{
  delete[] data;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlockIter member functions

SCMatrixBlockIter::SCMatrixBlockIter()
{
}

SCMatrixBlockIter::~SCMatrixBlockIter()
{
}

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
  void* casts[] =  { SavableState::_castdown(cd) };
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

/////////////////////////////////////////////////////////////////////////////
// SCDimension reference member functions

SavableState_REF_def(SCDimension);

/////////////////////////////////////////////////////////////////////////////
// SCMatrix reference base class member functions

SavableState_named_REF_def(RefSSSCMatrix,SCMatrix);
SavableState_named_REF_def(RefSSSymmSCMatrix,SymmSCMatrix);
SavableState_named_REF_def(RefSSDiagSCMatrix,DiagSCMatrix);

/////////////////////////////////////////////////////////////////////////////
// SCMatrix reference member functions
RefSCMatrix::RefSCMatrix() {}
             
RefSCMatrix::RefSCMatrix (RefSCMatrix & o): RefSSSCMatrix (o) {}
             
RefSCMatrix::RefSCMatrix (SCMatrix * o): RefSSSCMatrix (o) {}
             
RefSCMatrix::RefSCMatrix (RefDescribedClassBase&o): RefSSSCMatrix (o) {}

RefSCMatrix::~RefSCMatrix () {}

RefSCMatrix&
RefSCMatrix::operator=(SCMatrix* cr)
{
  RefSSSCMatrix::operator=(cr);
  return *this;
}

RefSCMatrix&
RefSCMatrix::operator=( RefDescribedClassBase & c)
{
  RefSSSCMatrix::operator=(c);
  return *this;
}

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
RefSCMatrix::operator+(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix ret(rowdim(),coldim());
  
  ret->copy(pointer());
  ret->accumulate(a);

  return ret;
}

RefSCMatrix
RefSCMatrix::operator-(RefSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSCMatrix ret(rowdim(),coldim());
  
  ret->copy(a.pointer());
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
  ret->copy(pointer());
  ret->transpose_this();
  return ret;
}

RefSCMatrix
RefSCMatrix::i()
{
  require_nonnull();
  
  RefSCMatrix ret;
  ret = clone();
  ret->copy(pointer());
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

void
RefSCMatrix::copy(RefSCMatrix&a)
{
  require_nonnull();
  pointer()->copy(a.pointer());
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

///////////////////////////////////////////////////////////////////
// RefSymmSCMatrix members

RefSymmSCMatrix::RefSymmSCMatrix()
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
             
RefSymmSCMatrix::RefSymmSCMatrix (RefDescribedClassBase&o):
  RefSSSymmSCMatrix (o)
{
}

RefSymmSCMatrix::~RefSymmSCMatrix ()
{
}

RefSymmSCMatrix&
RefSymmSCMatrix::operator=(SymmSCMatrix* cr)
{
  RefSSSymmSCMatrix::operator=(cr);
  return *this;
}

RefSymmSCMatrix&
RefSymmSCMatrix::operator=( RefDescribedClassBase & c)
{
  RefSSSymmSCMatrix::operator=(c);
  return *this;
}

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

RefSymmSCMatrix
RefSymmSCMatrix::operator+(RefSymmSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSymmSCMatrix ret(dim());
  
  ret->copy(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefSymmSCMatrix
RefSymmSCMatrix::operator-(RefSymmSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefSymmSCMatrix ret(dim());
  
  ret->copy(a.pointer());
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
  ret->copy(pointer());
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

SCMatrixdouble
RefSymmSCMatrix::operator()(int i,int j)
{
  return SCMatrixdouble(pointer(),i,j);
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

void
RefSymmSCMatrix::copy(RefSymmSCMatrix&a)
{
  require_nonnull();
  pointer()->copy(a.pointer());
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

///////////////////////////////////////////////////////////////////
// RefDiagSCMatrix members

RefDiagSCMatrix::RefDiagSCMatrix()
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
             
RefDiagSCMatrix::RefDiagSCMatrix (RefDescribedClassBase&o):
  RefSSDiagSCMatrix (o)
{
}

RefDiagSCMatrix::~RefDiagSCMatrix ()
{
}

RefDiagSCMatrix&
RefDiagSCMatrix::operator=(DiagSCMatrix* cr)
{
  RefSSDiagSCMatrix::operator=(cr);
  return *this;
}

RefDiagSCMatrix&
RefDiagSCMatrix::operator=( RefDescribedClassBase & c)
{
  RefSSDiagSCMatrix::operator=(c);
  return *this;
}

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
  pointer()->set_element(i,i,a);
}

double
RefDiagSCMatrix::get_element(int i)
{
  require_nonnull();
  return pointer()->get_element(i,i);
}

RefDiagSCMatrix
RefDiagSCMatrix::operator+(RefDiagSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix ret(dim());
  
  ret->copy(pointer());
  ret->accumulate(a.pointer());

  return ret;
}

RefDiagSCMatrix
RefDiagSCMatrix::operator-(RefDiagSCMatrix&a)
{
  require_nonnull();
  a.require_nonnull();

  RefDiagSCMatrix ret(dim());
  
  ret->copy(a.pointer());
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
  ret->copy(pointer());
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

SCMatrixdouble
RefDiagSCMatrix::operator()(int i)
{
  return SCMatrixdouble(pointer(),i,i);
}

SCMatrixdouble
RefDiagSCMatrix::operator()(int i,int j)
{
  return SCMatrixdouble(pointer(),i,j);
}

RefDiagSCMatrix
RefDiagSCMatrix::clone()
{
  RefDiagSCMatrix r = dim()->create_diagmatrix();
  return r;
}

void
RefDiagSCMatrix::copy(RefDiagSCMatrix&a)
{
  require_nonnull();
  pointer()->copy(a.pointer());
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
