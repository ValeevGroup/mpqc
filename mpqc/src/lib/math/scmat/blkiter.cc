
#ifdef __GNUC__
#pragma implementation
#endif

#include <iostream.h>
#include <math/scmat/blkiter.h>
#include <math/scmat/block.h>

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlockIter member functions

SCMatrixBlockIter::~SCMatrixBlockIter()
{
}

void
SCMatrixBlockIter::accum(double a)
{
  set(get()+a);
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixRectBlockIter member functions

SCMatrixRectBlockIter::SCMatrixRectBlockIter(SCMatrixRectBlock*a):
  block(a)
{
  reset();
}

void
SCMatrixRectBlockIter::reset()
{
  block_index = 0;
  i_ = block->istart;
  j_ = block->jstart;
}

SCMatrixRectBlockIter::~SCMatrixRectBlockIter()
{
}

int
SCMatrixRectBlockIter::i()
{
  return i_;
}

int
SCMatrixRectBlockIter::j()
{
  return j_;
}

double
SCMatrixRectBlockIter::get()
{
  return block->data[block_index];
}

void
SCMatrixRectBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCMatrixRectBlockIter::operator int()
{
  return (i_ < block->iend && j_ < block->jend);
}

void
SCMatrixRectBlockIter::operator ++()
{
  j_++;
  if (j_ >= block->jend) {
      j_ = block->jstart;
      i_++;
    }
  block_index++;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixRectSubBlockIter member functions

SCMatrixRectSubBlockIter::SCMatrixRectSubBlockIter(SCMatrixRectSubBlock*a):
  block(a)
{
  reset();
}

void
SCMatrixRectSubBlockIter::reset()
{
  i_ = block->istart;
  j_ = block->jstart;
  block_index = i_ * block->istride + j_;
}

SCMatrixRectSubBlockIter::~SCMatrixRectSubBlockIter()
{
}

int
SCMatrixRectSubBlockIter::i()
{
  return i_;
}

int
SCMatrixRectSubBlockIter::j()
{
  return j_;
}

double
SCMatrixRectSubBlockIter::get()
{
  return block->data[block_index];
}

void
SCMatrixRectSubBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCMatrixRectSubBlockIter::operator int()
{
  return (i_ < block->iend && j_ < block->jend);
}

void
SCMatrixRectSubBlockIter::operator ++()
{
  j_++;
  block_index++;
  if (j_ >= block->jend) {
      j_ = block->jstart;
      i_++;
      block_index += block->istride - (block->jend - block->jstart);
    }
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixLTriBlockIter member functions

SCMatrixLTriBlockIter::SCMatrixLTriBlockIter(SCMatrixLTriBlock*a):
  block(a)
{
  reset();
}

void
SCMatrixLTriBlockIter::reset()
{
  block_index = 0;
  i_ = block->start;
  j_ = block->start;
}

SCMatrixLTriBlockIter::~SCMatrixLTriBlockIter()
{
}

int
SCMatrixLTriBlockIter::i()
{
  return i_;
}

int
SCMatrixLTriBlockIter::j()
{
  return j_;
}

double
SCMatrixLTriBlockIter::get()
{
  return block->data[block_index];
}

void
SCMatrixLTriBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCMatrixLTriBlockIter::operator int()
{
  return (i_ < block->end);
}

void
SCMatrixLTriBlockIter::operator ++()
{
  j_++;
  if (j_ > i_) {
      j_ = block->start;
      i_++;
    }
  block_index++;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixLTriSubBlockIter member functions

SCMatrixLTriSubBlockIter::SCMatrixLTriSubBlockIter(
    SCMatrixLTriSubBlock*a):
  block(a)
{
  reset();
}

void
SCMatrixLTriSubBlockIter::reset()
{
  i_ = block->istart;
  j_ = block->jstart;
  block_index = (i_*(i_+1)>>1) + j_;
}

SCMatrixLTriSubBlockIter::~SCMatrixLTriSubBlockIter()
{
}

int
SCMatrixLTriSubBlockIter::i()
{
  return i_;
}

int
SCMatrixLTriSubBlockIter::j()
{
  return j_;
}

double
SCMatrixLTriSubBlockIter::get()
{
  return block->data[block_index];
}

void
SCMatrixLTriSubBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCMatrixLTriSubBlockIter::operator int()
{
  return (i_ < block->iend);
}

void
SCMatrixLTriSubBlockIter::operator ++()
{
  j_++;
  block_index++;
  if (j_ > i_) {
      j_ = block->jstart;
      i_++;
      block_index += block->istart;
    }
  if (j_ >= block->jend) {
      j_ = block->jstart;
      block_index += i_;
      i_++;
    }
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixDiagBlockIter member functions

SCMatrixDiagBlockIter::SCMatrixDiagBlockIter(SCMatrixDiagBlock*a):
  block(a)
{
  reset();
}

void
SCMatrixDiagBlockIter::reset()
{
  block_index = 0;
  i_ = block->istart;
}

SCMatrixDiagBlockIter::~SCMatrixDiagBlockIter()
{
}

int
SCMatrixDiagBlockIter::i()
{
  return i_;
}

int
SCMatrixDiagBlockIter::j()
{
  return i_ + block->jstart - block->istart;
}

double
SCMatrixDiagBlockIter::get()
{
  return block->data[block_index];
}

void
SCMatrixDiagBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCMatrixDiagBlockIter::operator int()
{
  return (i_ < block->iend);
}

void
SCMatrixDiagBlockIter::operator ++()
{
  i_++;
  block_index++;
}

/////////////////////////////////////////////////////////////////////////////
// SCMatrixDiagSubBlockIter member functions

SCMatrixDiagSubBlockIter::SCMatrixDiagSubBlockIter(SCMatrixDiagSubBlock*a):
  block(a)
{
  reset();
}

void
SCMatrixDiagSubBlockIter::reset()
{
  block_index = block->offset;
  i_ = block->istart;
}

SCMatrixDiagSubBlockIter::~SCMatrixDiagSubBlockIter()
{
}

int
SCMatrixDiagSubBlockIter::i()
{
  return i_;
}

int
SCMatrixDiagSubBlockIter::j()
{
  return i_ + block->jstart - block->istart;
}

double
SCMatrixDiagSubBlockIter::get()
{
  return block->data[block_index];
}

void
SCMatrixDiagSubBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCMatrixDiagSubBlockIter::operator int()
{
  return (i_ < block->iend);
}

void
SCMatrixDiagSubBlockIter::operator ++()
{
  i_++;
  block_index++;
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorSimpleBlockIter member functions

SCVectorSimpleBlockIter::SCVectorSimpleBlockIter(SCVectorSimpleBlock*a):
  block(a)
{
  reset();
}

void
SCVectorSimpleBlockIter::reset()
{
  block_index = 0;
  i_ = block->istart;
}

SCVectorSimpleBlockIter::~SCVectorSimpleBlockIter()
{
}

int
SCVectorSimpleBlockIter::i()
{
  return i_;
}

int
SCVectorSimpleBlockIter::j()
{
  fprintf(stderr,"SCVectorSimpleBlockIter::j() attempted to find j value\n");
  abort();
  return 0;
}

double
SCVectorSimpleBlockIter::get()
{
  return block->data[block_index];
}

void
SCVectorSimpleBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCVectorSimpleBlockIter::operator int()
{
  return (i_ < block->iend);
}

void
SCVectorSimpleBlockIter::operator ++()
{
  i_++;
  block_index++;
}

/////////////////////////////////////////////////////////////////////////////
// SCVectorSimpleSubBlockIter member functions

SCVectorSimpleSubBlockIter::SCVectorSimpleSubBlockIter(
    SCVectorSimpleSubBlock*a):
  block(a)
{
  reset();
}

void
SCVectorSimpleSubBlockIter::reset()
{
  block_index = block->offset;
  i_ = block->istart;
}

SCVectorSimpleSubBlockIter::~SCVectorSimpleSubBlockIter()
{
}

int
SCVectorSimpleSubBlockIter::i()
{
  return i_;
}

int
SCVectorSimpleSubBlockIter::j()
{
  fprintf(stderr,"SCVectorSimpleSubBlockIter::j():"
          "attempted to find j value\n");
  abort();
  return 0;
}

double
SCVectorSimpleSubBlockIter::get()
{
  return block->data[block_index];
}

void
SCVectorSimpleSubBlockIter::set(double a)
{
  block->data[block_index] = a;
}

SCVectorSimpleSubBlockIter::operator int()
{
  return (i_ < block->iend);
}

void
SCVectorSimpleSubBlockIter::operator ++()
{
  i_++;
  block_index++;
}
