//
// tensor.cc
//
// Copyright (C) 2008 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@qtp.ufl.edu>
// Maintainer: TS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdexcept>
#include <sstream>
#include <cassert>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>
#include <util/misc/string.h>
#include <util/class/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>
//F77 blas routines
#include <chemistry/qc/mbptr12/blas.h>

#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

using std::map;
using std::pair;
using std::vector;
using namespace std;
using namespace sc;

static ClassDesc Tensor_cd(
  typeid(Tensor), "Tensor", 1, "virtual public RefCount",
  0, 0, 0);

Tensor::Tensor(std::string filename, const Ref<MemoryGrp>& mem): mem_(mem) {
  filename_ = ".smith_" + filename;
  file_allocated_ = false;
}

Tensor::~Tensor(){
  assert(file_allocated_);
  deletefile();
}


void Tensor::input_offset(long tag, long offset){
  pair<std::map<long, long>::iterator, bool> val = hash_table_.insert(make_pair(tag, offset));
  assert(val.second);
}


void Tensor::set_filesize(long i){
  assert(hash_table_.size() > 0);
  filesize_ = i;
  long tag = LONG_MAX / 100; // infty; divided by 100 so that it will not raise any problems
  hash_table_.insert(std::map<long, long>::value_type(tag, i));
}


void Tensor::createfile(){
 mem_->sync();
 vector<long> filesizes = determine_filesizes();
 long maxsize = *max_element(filesizes.begin(), filesizes.end());

 file_ = new MemoryGrpRegion(mem_, sizeof(double) * maxsize);  // should be modified
 file()->set_localsize(sizeof(double) * filesizes[mem_->me()]);

 file_allocated_ = true;
 zero();
}


vector<long> Tensor::determine_filesizes(){
 vector<long> filesizes;
 long defaultsize = filesize_ / mem_->n();
 if (filesize_ % mem_->n() != 0) ++defaultsize;

 long current;
 long previous = 0L;
 for (long i = 1L; i < mem_->n(); ++i){
  long thresh = min(defaultsize * i, filesize_);

  std::map<long, long>::iterator j = hash_table_.begin();
  while (j->second < thresh) ++j;

  if (j == hash_table_.begin()) {
    std::map<long, long>::iterator k = j; k++;
    current = thresh < (k->second) - thresh ? hash_table_.begin()->second : k->second;
  } else if (j == (--hash_table_.end())) {
    current = j->second;
    assert(current == filesize_);
  } else {
    std::map<long, long>::iterator k = j; k++;
    current = thresh - (j->second) <= (k->second) - thresh ? j->second : k->second;
  }
  filesizes.push_back(current - previous);
  previous = current;
 }
 filesizes.push_back(filesize_ - previous);
 return filesizes;
}


void Tensor::deletefile(){
 delete file_;
}


/// routines called from one node (i.e. inside the loops) >>>>>>>>>>>>>>>>>>>>>>>>>
void Tensor::get_block(long tag, double* data){
  std::map<long, long>::iterator iter = hash_table_.find(tag);
  assert(iter != hash_table_.end());
  long doffset = iter->second;
  long dsize   = (++iter)->second-doffset;
  distsize_t offset = (distsize_t)doffset * sizeof(double);
  const int    size =        (int)dsize   * sizeof(double); // in byte

  double* buffer = (double*) file()->obtain_readonly(offset, size);
  ::memcpy((void*) data, (void*) buffer, size);
  file()->release_readonly(buffer, offset, size);
}


void Tensor::put_block(long tag, double* data){
  std::map<long,long>::iterator iter = hash_table_.find(tag);
  assert(iter != hash_table_.end());
  long doffset = iter->second;
  long dsize   = (++iter)->second - doffset;
  distsize_t offset = (distsize_t)doffset * sizeof(double);
  const int    size =        (int)dsize   * sizeof(double); // in byte

  double* buffer = (double *) file()->obtain_writeonly(offset, size);
  ::memcpy((void*) buffer, (void*) data, size);
  file()->release_writeonly((void*) buffer, offset, size);
}


void Tensor::add_block(long tag, double* data){
// copy of data will be allocated;
// probably it is ok (in smith codes, k_a and k_a_sort are
// already deallocated at the time this is called)...
  std::map<long,long>::iterator iter = hash_table_.find(tag);
  assert(iter != hash_table_.end());
  long doffset = iter->second;
  int dsize    = (int)((++iter)->second - doffset);
  distsize_t offset = (distsize_t) doffset * sizeof(double);
  const int    size =          (int)dsize  * sizeof(double); // in byte

  double* copy_data = mem_->malloc_local_double(dsize);
  ::memcpy((void*) copy_data, (void*) data, size);

  double* buffer = (double*) file()->obtain_readonly(offset, size);
  const int unit = 1;
  const double one = 1.0;
  F77_DAXPY(&dsize, &one, buffer, &unit, copy_data, &unit);
  file()->release_readonly((void*) buffer, offset, size);

  buffer = (double*) file()->obtain_writeonly(offset, size);
  ::memcpy((void*) buffer, (void*) copy_data, size);
  file()->release_writeonly((void*) buffer, offset, size);

  mem_->free_local_double(copy_data);
}

bool Tensor::exists(long tag) const {
  std::map<long, long>::const_iterator iter = hash_table_.find(tag);
  if (iter == hash_table_.end()) return false;
  return true;
}

bool Tensor::is_this_local(long tag){
  std::map<long,long>::iterator iter = hash_table_.find(tag);
  if (iter == hash_table_.end()) return false;

  long myoffset_   = iter->second * sizeof(double); // in byte
  long localoffset_= (long)(file_->localoffset());
  long nextoffset_ = (long)(file_->localsize()) + localoffset_;
  return localoffset_ <= myoffset_ && myoffset_ < nextoffset_;
}


/// routines called from all the nodes at the same time >>>>>>>>>>>>>>>>>>>>>>>>>
void Tensor::zero(){
  assign(0.0);
}

void Tensor::assign(double a){
  double* buffer = (double *) file()->localdata();
  const size_t dsize = (size_t)file()->localsize() / sizeof(double);
  std::fill(buffer, buffer + dsize, a);
  sync();
}

namespace {
  // UnaryFunction scales its argument by a constant
  struct ElementScaler : public std::unary_function<double&, void> {
    public:
      ElementScaler(double scale) : scale_(scale) {}
      void operator()(double& v) { v *= scale_; }
    private:
      double scale_;
  };
}
void Tensor::scale(double a){
  double* buffer=(double *) file()->localdata();
  const size_t dsize =(size_t)file()->localsize() / sizeof(double);
  ElementScaler scaler(a);
  std::for_each(buffer, buffer + dsize, scaler);
  sync();
}

void Tensor::daxpy(const Ref<Tensor>& other, double a){ // add to self
  const int dsize = file()->localsize() / sizeof(double);
  assert(dsize == other->file()->localsize() / sizeof(double));
  double* buffer1 = (double *) file()->localdata();
  const double* buffer2 = (double *) other->file()->localdata();
  const int unit = 1;
  F77_DAXPY(&dsize, &a, buffer2, &unit, buffer1, &unit);
  sync();
}


Ref<Tensor> Tensor::copy() const{
// TODO this is not sufficient for the disk-based algorithm
// perhaps we need to check filenames in the current directory, and find
// an appropriate name for the tensor to be created..
  Ref<Tensor> other = new Tensor("TODO", mem_);
  other->hash_table_ = hash_table_;
  other->filesize_   = filesize_;
  other->createfile();

  int size = file()->localsize();
  assert(size == other->file()->localsize());
  double* buffer1 = (double *) file()->localdata();
  double* buffer2 = (double *) other->file()->localdata();
  ::memcpy((void*) buffer2, (void*) buffer1, size);
  sync();

  return other;
}


Ref<Tensor> Tensor::clone() const {
  Ref<Tensor> other = new Tensor("TODO",mem_);
  other->hash_table_ = hash_table_;
  other->filesize_   = filesize_;
  other->createfile();
  return other;
}


double Tensor::norm() const {
  sync();
  const double* buffer = (double *) file()->localdata();
  const int dsize = file()->localsize() / sizeof(double);
  const int unit = 1;
  double norm_ = F77_DDOT(&dsize, buffer, &unit, buffer, &unit);

  Ref<MessageGrp> msg_ = MessageGrp::get_default_messagegrp();
  msg_->sum(norm_);

  /// sqrt needed?
  norm_ = ::sqrt(norm_);
  return norm_;
}


double Tensor::ddot(Ref<Tensor>& other) const {
  sync();
  const double* buffer1 = (double *)        file()->localdata();
  const double* buffer2 = (double *) other->file()->localdata();
  const int dsize = file()->localsize() / sizeof(double);
  assert(file()->localsize() == other->file()->localsize());
  const int unit = 1;
  double ddotproduct = F77_DDOT(&dsize, buffer1, &unit, buffer2, &unit);

  sync();
  Ref<MessageGrp> msg_ = MessageGrp::get_default_messagegrp();
  msg_->sum(ddotproduct);
  return ddotproduct;
}

double sc::RMS(const Tensor& t) {
  return t.norm() / t.get_filesize();
}

