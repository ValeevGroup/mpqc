
#include <math.h>
#include <math/scmat/abstract.h>
#include <math/scmat/blkiter.h>

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCMatrix classes.
/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
// SCDimension members

#define CLASSNAME SCDimension
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCDimension::SCDimension()
{
}

void *
SCDimension::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

SCDimension::~SCDimension()
{
}

/////////////////////////////////////////////////////////////////////////
// SCElementScale members

#define CLASSNAME SCElementScale
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementScale::SCElementScale(double a):scale(a) {}
SCElementScale::SCElementScale(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(scale);
}
void
SCElementScale::save_data_state(StateOut&s)
{
  s.put(scale);
}
void *
SCElementScale::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
  return do_castdowns(casts,cd);
}
SCElementScale::~SCElementScale() {}
void
SCElementScale::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(scale*i.get());
    }
}

/////////////////////////////////////////////////////////////////////////
// SCElementMaxAbs members

SavableState_REF_def(SCElementMaxAbs);
#define CLASSNAME SCElementMaxAbs
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

SCElementMaxAbs::SCElementMaxAbs(double a):r(0.0) {}
SCElementMaxAbs::SCElementMaxAbs(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(r);
}
void
SCElementMaxAbs::save_data_state(StateOut&s)
{
  s.put(r);
}
void *
SCElementMaxAbs::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
  return do_castdowns(casts,cd);
}
SCElementMaxAbs::~SCElementMaxAbs() {}
void
SCElementMaxAbs::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (fabs(i.get()) > r) r = fabs(i.get());
    }
}
double
SCElementMaxAbs::result()
{
  return r;
}
int
SCElementMaxAbs::has_collect()
{
  return 1;
}
void
SCElementMaxAbs::collect(RefSCElementOp&op)
{
  RefSCElementMaxAbs ma(op);
  if (ma->r > r) r = ma->r;
}

/////////////////////////////////////////////////////////////////////////
// SCElementAssign members

#define CLASSNAME SCElementAssign
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementAssign::SCElementAssign(double a):assign(a) {}
SCElementAssign::SCElementAssign(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(assign);
}
void
SCElementAssign::save_data_state(StateOut&s)
{
  s.put(assign);
}
void *
SCElementAssign::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
  return do_castdowns(casts,cd);
}
SCElementAssign::~SCElementAssign() {}
void
SCElementAssign::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      i.set(assign);
    }
}

/////////////////////////////////////////////////////////////////////////
// SCElementShiftDiagonal members

#define CLASSNAME SCElementShiftDiagonal
#define PARENTS       virtual public SCDiagElementOp, \
                      virtual public SCSymmElementOp, \
                      virtual public SCRectElementOp, \
                      virtual public SCVectorElementOp
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
SCElementShiftDiagonal::SCElementShiftDiagonal(double a):shift_diagonal(a) {}
SCElementShiftDiagonal::SCElementShiftDiagonal(StateIn&s):
  SavableState(s,class_desc_)
{
  s.get(shift_diagonal);
}
void
SCElementShiftDiagonal::save_data_state(StateOut&s)
{
  s.put(shift_diagonal);
}
void *
SCElementShiftDiagonal::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SCDiagElementOp::_castdown(cd),
                     SCSymmElementOp::_castdown(cd),
                     SCRectElementOp::_castdown(cd),
                     SCVectorElementOp::_castdown(cd)
                   };
  return do_castdowns(casts,cd);
}
SCElementShiftDiagonal::~SCElementShiftDiagonal() {}
void
SCElementShiftDiagonal::process(SCMatrixBlockIter&i)
{
  for (i.reset(); i; ++i) {
      if (i.i() == i.j()) i.set(shift_diagonal+i.get());
    }
}

/////////////////////////////////////////////////////////////////////////
// SCMatrix members

#define CLASSNAME SCMatrix
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCMatrix::SCMatrix()
{
}

SCMatrix::~SCMatrix()
{
}

SCMatrix::SCMatrix(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
SCMatrix::save_data_state(StateOut&s)
{
}

void *
SCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
SCMatrix::nrow()
{
  return rowdim()->n();
}

int
SCMatrix::ncol()
{
  return coldim()->n();
}

double
SCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCRectElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SCMatrix::assign(double a)
{
  RefSCRectElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCMatrix::scale(double a)
{
  RefSCRectElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SCMatrix::shift_diagonal(double a)
{
  RefSCRectElementOp op = new SCElementShiftDiagonal(a);
  this->element_op(op);
}

void
SCMatrix::unit()
{
  this->assign(0.0);
  this->shift_diagonal(1.0);
}

void
SCMatrix::copy(SCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SCMatrix::print(ostream&o)
{
  print(0, o, 10);
}

/////////////////////////////////////////////////////////////////////////
// SymmSCMatrix member functions

#define CLASSNAME SymmSCMatrix
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SymmSCMatrix::SymmSCMatrix()
{
}

SymmSCMatrix::SymmSCMatrix(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
SymmSCMatrix::save_data_state(StateOut&s)
{
}

void *
SymmSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
SymmSCMatrix::n()
{
  return dim()->n();
}

double
SymmSCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCSymmElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SymmSCMatrix::assign(double a)
{
  RefSCSymmElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SymmSCMatrix::scale(double a)
{
  RefSCSymmElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SymmSCMatrix::shift_diagonal(double a)
{
  RefSCSymmElementOp op = new SCElementShiftDiagonal(a);
  this->element_op(op);
}

void
SymmSCMatrix::unit()
{
  this->assign(0.0);
  this->shift_diagonal(1.0);
}

void
SymmSCMatrix::copy(SymmSCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SymmSCMatrix::print(ostream&o)
{
  print(0, o, 10);
}

/////////////////////////////////////////////////////////////////////////
// DiagSCMatrix member functions

#define CLASSNAME DiagSCMatrix
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

DiagSCMatrix::DiagSCMatrix()
{
}

DiagSCMatrix::DiagSCMatrix(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
DiagSCMatrix::save_data_state(StateOut&s)
{
}

void *
DiagSCMatrix::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
DiagSCMatrix::n()
{
  return dim()->n();
}

double
DiagSCMatrix::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCDiagElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
DiagSCMatrix::assign(double a)
{
  RefSCDiagElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
DiagSCMatrix::scale(double a)
{
  RefSCDiagElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
DiagSCMatrix::copy(DiagSCMatrix*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
DiagSCMatrix::print(ostream&o)
{
  print(0, o, 10);
}

/////////////////////////////////////////////////////////////////////////
// These member are used by the abstract SCVector classes.
/////////////////////////////////////////////////////////////////////////

#define CLASSNAME SCVector
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SCVector::SCVector()
{
}

SCVector::~SCVector()
{
}

SCVector::SCVector(StateIn&s):
  SavableState(s,class_desc_)
{
}

void
SCVector::save_data_state(StateOut&s)
{
}

void *
SCVector::_castdown(const ClassDesc*cd)
{
  void* casts[] =  { SavableState::_castdown(cd) };
  return do_castdowns(casts,cd);
}

int
SCVector::n()
{
  return dim()->n();
}

double
SCVector::maxabs()
{
  RefSCElementMaxAbs op = new SCElementMaxAbs();
  RefSCRectElementOp abop = op;
  this->element_op(abop);
  return op->result();
}

void
SCVector::assign(double a)
{
  RefSCVectorElementOp op = new SCElementAssign(a);
  this->element_op(op);
}

void
SCVector::scale(double a)
{
  RefSCVectorElementOp op = new SCElementScale(a);
  this->element_op(op);
}

void
SCVector::copy(SCVector*a)
{
  this->assign(0.0);
  this->accumulate(a);
}

void
SCVector::print(ostream&o)
{
  print(0, o, 10);
}
