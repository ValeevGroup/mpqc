
#ifndef _util_group_memproc_cc
#define _util_group_memproc_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/memproc.h>

#define CLASSNAME ProcMemoryGrp
#define PARENTS public MemoryGrp
#include <util/class/classi.h>
void *
ProcMemoryGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MemoryGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ProcMemoryGrp::ProcMemoryGrp(int localsize)
{
  offsets_ = new int[2];
  offsets_[0] = 0;
  offsets_[1] = localsize;
  n_ = 1;
  me_ = 0;
  data_ = new char[localsize];
}

ProcMemoryGrp::~ProcMemoryGrp()
{
  delete[] data_;
}

void *
ProcMemoryGrp::obtain_readwrite(int offset, int size)
{
  return &data_[offset];
}

void *
ProcMemoryGrp::obtain_readonly(int offset, int size)
{
  return &data_[offset];
}

void
ProcMemoryGrp::release_read(void *data, int offset, int size)
{
}

void
ProcMemoryGrp::release_write(void *data, int offset, int size)
{
}

void
ProcMemoryGrp::sync()
{
}

#endif
