
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/blocked.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/elemop.h>

/////////////////////////////////////////////////////////////////////////////
// BlockedSCMatrixKit member functions

SavableState_REF_def(BlockedSCMatrixKit);

#define CLASSNAME BlockedSCMatrixKit
#define PARENTS public SCMatrixKit
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
BlockedSCMatrixKit::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCMatrixKit::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCMatrixKit::BlockedSCMatrixKit(const RefSCMatrixKit&subkit):
  subkit_(subkit)
{
}

BlockedSCMatrixKit::BlockedSCMatrixKit(const RefKeyVal& keyval):
  SCMatrixKit(keyval)
{
  subkit_ = keyval->describedclassvalue("subkit");
}

BlockedSCMatrixKit::~BlockedSCMatrixKit()
{
}

SCMatrix*
BlockedSCMatrixKit::matrix(const RefSCDimension&d1, const RefSCDimension&d2)
{
  int i;
  for (i=0; i<d1->blocks()->nblock(); i++) {
      if (d1->blocks()->subdim(i).null()) {
          cerr << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  for (i=0; i<d2->blocks()->nblock(); i++) {
      if (d2->blocks()->subdim(i).null()) {
          cerr << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedSCMatrix(d1,d2,this);
}

SymmSCMatrix*
BlockedSCMatrixKit::symmmatrix(const RefSCDimension&d)
{
  for (int i=0; i<d->blocks()->nblock(); i++) {
      if (d->blocks()->subdim(i).null()) {
          cerr << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedSymmSCMatrix(d,this);
}

DiagSCMatrix*
BlockedSCMatrixKit::diagmatrix(const RefSCDimension&d)
{
  for (int i=0; i<d->blocks()->nblock(); i++) {
      if (d->blocks()->subdim(i).null()) {
          cerr << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedDiagSCMatrix(d,this);
}

SCVector*
BlockedSCMatrixKit::vector(const RefSCDimension&d)
{
  for (int i=0; i<d->blocks()->nblock(); i++) {
      if (d->blocks()->subdim(i).null()) {
          cerr << indent
               << "BlockedSCMatrixKit: given a dim without subdim info"
               << endl;
          abort();
        }
    }
  return new BlockedSCVector(d,this);
}

/////////////////////////////////////////////////////////////////////////////

#define CLASSNAME BlockedSCElementOp
#define PARENTS public SCElementOP
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
BlockedSCElementOp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCElementOp::BlockedSCElementOp()
{
  current_block_=0;
}

void
BlockedSCElementOp::working_on(int b)
{
  current_block_ = b;
}

int
BlockedSCElementOp::current_block() const
{
  return current_block_;
}

#define CLASSNAME BlockedSCElementOp2
#define PARENTS public SCElementOP2
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
BlockedSCElementOp2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp2::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCElementOp2::BlockedSCElementOp2()
{
  current_block_=0;
}

void
BlockedSCElementOp2::working_on(int b)
{
  current_block_ = b;
}

int
BlockedSCElementOp2::current_block() const
{
  return current_block_;
}

#define CLASSNAME BlockedSCElementOp3
#define PARENTS public SCElementOP3
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
BlockedSCElementOp3::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SCElementOp3::_castdown(cd);
  return do_castdowns(casts,cd);
}

BlockedSCElementOp3::BlockedSCElementOp3()
{
  current_block_=0;
}

void
BlockedSCElementOp3::working_on(int b)
{
  current_block_ = b;
}

int
BlockedSCElementOp3::current_block() const
{
  return current_block_;
}
