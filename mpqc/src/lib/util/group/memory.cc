
#ifndef _util_group_memory_cc
#define _util_group_memory_cc

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/group/memory.h>

//////////////////////////////////////////////////////////////////////
// MemoryGrpBuf members

template <class data_t>
MemoryGrpBuf<data_t>::MemoryGrpBuf(const RefMemoryGrp & grp)
{
  grp_ = grp;
  locktype_ = None;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::writeonly(int offset, int length)
{
  if (locktype_ != None) release();
  data_ = (data_t *) grp_->obtain_writeonly(sizeof(data_t)*offset,
                                            sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  locktype_ = Write;
  return data_;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::readwrite(int offset, int length)
{
  if (locktype_ != None) release();
  data_ = (data_t *) grp_->obtain_readwrite(sizeof(data_t)*offset,
                                            sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  locktype_ = Write;
  return data_;
}

template <class data_t>
const data_t *
MemoryGrpBuf<data_t>::readonly(int offset, int length)
{
  if (locktype_ != None) release();
  data_ = (data_t *) grp_->obtain_readonly(sizeof(data_t)*offset,
                                           sizeof(data_t)*length);
  offset_ = offset;
  length_ = length;
  locktype_ = Read;
  return data_;
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::writeonly_on_node(int offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return writeonly(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
data_t *
MemoryGrpBuf<data_t>::readwrite_on_node(int offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return readwrite(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
const data_t *
MemoryGrpBuf<data_t>::readonly_on_node(int offset, int length, int node)
{
  if (node == -1) node = grp_->me();
  return readonly(offset + grp_->offset(node)/sizeof(data_t), length);
}

template <class data_t>
void
MemoryGrpBuf<data_t>::release()
{
  if (locktype_ == Write)
      grp_->release_write((data_t *)data_,
                          sizeof(data_t)*offset_, sizeof(data_t)*length_);
  if (locktype_ == Read)
      grp_->release_read(data_, sizeof(data_t)*offset_, sizeof(data_t)*length_);

  locktype_ = None;
}

#endif
