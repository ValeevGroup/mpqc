
#include <math/scmat/blkiter.h>
#include <math/scmat/block.h>

/////////////////////////////////////////////////////////////////////////////
// SCMatrixBlockIter member functions

SCMatrixBlockIter::SCMatrixBlockIter()
{
}

SCMatrixBlockIter::~SCMatrixBlockIter()
{
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
