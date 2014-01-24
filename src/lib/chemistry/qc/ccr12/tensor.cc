//
// tensor.cc
//
// Copyright (C) 2008 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

#include <climits>
#include <cassert>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <cstdio>
#include <unistd.h>
#include <util/misc/exenv.h>
#include <util/misc/exenv.h>
#include <math/scmat/blas.h>
#include <util/misc/consumableresources.h>
#include <chemistry/qc/ccr12/tensor.h>
#include <chemistry/qc/ccr12/ccr12_info.h>

using namespace std;
using namespace sc;

static ClassDesc Tensor_cd(
  typeid(Tensor), "Tensor", 1, "virtual public RefCount",
  0, 0, 0);

Tensor::Tensor(string filename, const Ref<MemoryGrp>& mem): mem_(mem) {
  const std::string basename_prefix = SCFormIO::fileext_to_filename(".smith.");
  const std::string full_prefix = ConsumableResources::get_default_instance()->disk_location() +
                                  basename_prefix;
  filename_ = full_prefix + filename;
  file_allocated_ = false;
}

Tensor::~Tensor(){
  MPQC_ASSERT(file_allocated_);
  deletefile();
}


void Tensor::input_offset(long tag, long offset){
  hash_table_.insert(std::map<long, long>::value_type(tag, offset));
}


void Tensor::set_filesize(long i){
  MPQC_ASSERT(hash_table_.size() > 0);
  filesize_ = i;
  long tag = LONG_MAX - 1;
  hash_table_.insert(std::map<long, long>::value_type(tag, i));
}


#ifdef DISK_BASED_SMITH
// must fit into the cache of the hard drive (?)
static const long cachesize = 100000;
#endif

void Tensor::createfile(){
#ifndef DISK_BASED_SMITH
 mem_->sync();
 vector<long> filesizes = determine_filesizes();
 long maxsize = *max_element(filesizes.begin(), filesizes.end());

 file_ = new MemoryGrpRegion(mem_, sizeof(double) * maxsize);
 file_->set_localsize(sizeof(double) * filesizes[mem_->me()]);

 zero();
#else
 // create a zero-cleared file of the specified size
 file_ = new fstream(filename_.c_str(), fstream::out | fstream::binary); 
 MPQC_ASSERT(file()->is_open());
 double* aux_array = new double[cachesize];
 fill(aux_array, aux_array + cachesize, 0.0);
 long size_back = filesize_; 
 while (size_back > 0L) {
   file_->write((const char*)aux_array, min(cachesize, size_back) * sizeof(double)); 
   size_back -= cachesize;
 }
 delete[] aux_array;
 delete file_;

 // reopen with in/out/binary; will be delete'd in the destructor
 file_ = new fstream(filename_.c_str(), fstream::in | fstream::out | fstream::binary); 
 // Unlinking while open to avoid garbage. When closed, the file is physically removed.
 unlink(filename_.c_str());
#endif

 file_allocated_ = true;
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
    MPQC_ASSERT(current == filesize_);
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
  MPQC_ASSERT(iter != hash_table_.end());
  long doffset = iter->second;
  long dsize   = (++iter)->second-doffset;
  distsize_t offset = (distsize_t)doffset * sizeof(double);
  const int    size =        (int)dsize   * sizeof(double); // in byte

#ifndef DISK_BASED_SMITH
  double* buffer = (double*) file_->obtain_readonly(offset, size);
  ::memcpy((void*) data, (void*) buffer, size);
  file_->release_readonly(buffer, offset, size);
#else
  file_->clear();
  file_->seekg(doffset * sizeof(double));
  file_->read((char*) data, size);
#endif
}


void Tensor::put_block(long tag, double* data){
  std::map<long,long>::iterator iter = hash_table_.find(tag);
  MPQC_ASSERT(iter != hash_table_.end());
  long doffset = iter->second;
  long dsize   = (++iter)->second - doffset;
  distsize_t offset = (distsize_t)doffset * sizeof(double);
  const int    size =        (int)dsize   * sizeof(double); // in byte

#ifndef DISK_BASED_SMITH
  double* buffer = (double *) file_->obtain_writeonly(offset, size);
  ::memcpy((void*) buffer, (void*) data, size);
  file_->release_writeonly((void*) buffer, offset, size);
#else
  file_->clear();
  file_->seekp(doffset * sizeof(double));
  file_->write((char*)data, size);
#endif
}


void Tensor::add_block(long tag, double* data){
// copy of data will be allocated;
// probably it is ok (in smith codes, k_a and k_a_sort are
// already deallocated at the time this is called)...
  std::map<long,long>::iterator iter = hash_table_.find(tag);
  MPQC_ASSERT(iter != hash_table_.end());
  long doffset = iter->second;
  blasint dsize    = (int)((++iter)->second - doffset);
  distsize_t offset = (distsize_t) doffset * sizeof(double);
  const int    size =          (int)dsize  * sizeof(double); // in byte

  double* copy_data = mem_->malloc_local_double(dsize);
#ifndef DISK_BASED_SMITH
  ::memcpy((void*) copy_data, (void*) data, size);

  double* buffer = (double*) file_->obtain_readonly(offset, size);
  const blasint unit = 1;
  const double one = 1.0;
  F77_DAXPY(&dsize, &one, buffer, &unit, copy_data, &unit);
  file_->release_readonly((void*) buffer, offset, size);

  buffer = (double*) file_->obtain_writeonly(offset, size);
  ::memcpy((void*) buffer, (void*) copy_data, size);
  file_->release_writeonly((void*) buffer, offset, size);
#else
  file_->clear();
  file_->seekg(doffset * sizeof(double));
  file_->read((char*)copy_data, size);
  const double one = 1.0;
  const blasint unit = 1;
  F77_DAXPY(&dsize, &one, data, &unit, copy_data, &unit);

  file_->clear();
  file_->seekp(doffset * sizeof(double));
  file_->write((const char*)copy_data, size);
#endif

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

#ifndef DISK_BASED_SMITH
  long myoffset_   = iter->second * sizeof(double); // in byte
  long localoffset_= (long)(file()->localoffset());
  long nextoffset_ = (long)(file()->localsize()) + localoffset_;
  return localoffset_ <= myoffset_ && myoffset_ < nextoffset_;
#else
  // Only serial runs are supported with DISK_BASED_SMITH so far
  return true;
#endif
}


/// routines called from all the nodes at the same time >>>>>>>>>>>>>>>>>>>>>>>>>
void Tensor::zero(){
  assign(0.0);
}

void Tensor::assign(double a){
#ifndef DISK_BASED_SMITH
  double* buffer = (double *) file_->localdata();
  const size_t dsize = (size_t)file()->localsize() / sizeof(double);
  fill(buffer, buffer + dsize, a);
  sync();
#else
  file_->clear();
  file_->seekp(0);
  double* aux_array = new double[cachesize];
  fill(aux_array, aux_array + cachesize, a);
  long size_back = filesize_;
  while (size_back > 0) {
    const size_t readsize = min(size_back, cachesize) * sizeof(double);
    file_->write((const char*)aux_array, readsize); 
    size_back -= cachesize;
  }
  delete[] aux_array;
#endif
}

namespace {
  // UnaryFunction scales its argument by a constant
  struct ElementScaler : public unary_function<double&, void> {
    public:
      ElementScaler(double scale) : scale_(scale) {}
      void operator()(double& v) { v *= scale_; }
    private:
      double scale_;
  };
}
void Tensor::scale(double a){
#ifndef DISK_BASED_SMITH
  double* buffer=(double *) file_->localdata();
  const size_t dsize =(size_t)file()->localsize() / sizeof(double);
  ElementScaler scaler(a);
  for_each(buffer, buffer + dsize, scaler);
  sync();
#else
  file_->clear();
  double* aux_array = new double [cachesize];
  size_t size_now = 0UL; 
  long size_back = filesize_;
  ElementScaler scaler(a);
  while (size_back > 0) {
   const int rsize = min(cachesize, size_back);
   const size_t readsize = rsize * sizeof(double);
   const size_t position = size_now * sizeof(double);

   file_->seekg(position);
   file_->read((char*)aux_array, readsize); 
   for_each(aux_array, aux_array + rsize, scaler);

   file_->seekp(position);
   file_->write((const char*)aux_array, readsize); 

   size_now += cachesize;
   size_back -= cachesize;
 }
 delete[] aux_array;
#endif
}

void Tensor::daxpy(const Ref<Tensor>& other, double a){ // add to self
  const blasint unit = 1;
#ifndef DISK_BASED_SMITH
  const blasint dsize = file_->localsize() / sizeof(double);
  MPQC_ASSERT(dsize == other->file()->localsize() / sizeof(double));
  double* buffer1 = (double *) file_->localdata();
  const double* buffer2 = (double *) other->file()->localdata();
  F77_DAXPY(&dsize, &a, buffer2, &unit, buffer1, &unit);
  sync();
#else
// Assuming that this and other point to different Tensor's. Probably OK, right?
// Otherwise, we need to insert fstream::clear for the last block.
  double* aux_array = new double[cachesize];
  double* aux_array2 = new double[cachesize];
  size_t size_now = 0UL; 
  long size_back = filesize_;
  this->file_->clear();
  other->file_->clear();
  while (size_back > 0L) {
    const blasint rsize = min(cachesize, size_back);
    const size_t readsize = rsize * sizeof(double);
    const size_t position = size_now * sizeof(double); 

    this->file_->seekg(position);
    this->file_->read((char*)aux_array, readsize);
    other->file_->seekg(position);
    other->file_->read((char*)aux_array2, readsize);
    F77_DAXPY(&rsize, &a, aux_array2, &unit, aux_array, &unit);

    this->file_->seekp(position);
    this->file_->write((const char*)aux_array, readsize); 

    size_now += cachesize;
    size_back -= cachesize;
  }
  delete[] aux_array;
  delete[] aux_array2;
#endif
}


// This is an identifier for temp files
static int filecounter = 0;

Ref<Tensor> Tensor::copy() const{
  stringstream ss; 
  ss << "temp_file_" << filecounter;
  ++filecounter;
  Ref<Tensor> other = new Tensor(ss.str(), mem_);
  other->hash_table_ = hash_table_;
  other->filesize_   = filesize_;
  other->createfile();

#ifndef DISK_BASED_SMITH
  int size = file_->localsize();
  MPQC_ASSERT(size == other->file()->localsize());
  double* buffer1 = (double *) file_->localdata();
  double* buffer2 = (double *) other->file()->localdata();
  ::memcpy((void*) buffer2, (void*) buffer1, size);
  sync();
#else
  file_->clear();
  other->file()->clear();
  file_->seekg(0);
  other->file()->seekp(0);
  char* aux_array = new char[cachesize * sizeof(double)];
  long size_back = filesize_; 
  while (size_back > 0L) {
    const size_t readsize = min(cachesize, size_back) * sizeof(double);
    file_->read(aux_array, readsize);
    other->file()->write(aux_array, readsize);
    size_back -= cachesize;
  }
  delete[] aux_array;
#endif

  return other;
}


Ref<Tensor> Tensor::clone() const {
  stringstream ss; 
  ss << "temp_file_" << filecounter;
  ++filecounter;
  Ref<Tensor> other = new Tensor(ss.str(), mem_);
  other->hash_table_ = hash_table_;
  other->filesize_   = filesize_;
  other->createfile();
  return other;
}


double Tensor::norm() const {
  const blasint unit = 1;
#ifndef DISK_BASED_SMITH 
  sync();
  const double* buffer = (double *) file_->localdata();
  const blasint dsize = file_->localsize() / sizeof(double);
  double norm_ = F77_DDOT(&dsize, buffer, &unit, buffer, &unit);

  Ref<MessageGrp> msg_ = MessageGrp::get_default_messagegrp();
  msg_->sum(norm_);
#else
  file_->clear();
  file_->seekg(0);
  double* aux_array = new double[cachesize];
  double norm_ = 0.0;
  long size_back = filesize_;
  while (size_back > 0L) {
    const blasint bsize = min(cachesize, size_back);
    file_->read((char*)aux_array, bsize * sizeof(double));
    norm_ += F77_DDOT(&bsize, aux_array, &unit, aux_array, &unit);
    size_back -= cachesize;
  }
  delete[] aux_array;
#endif

  /// sqrt needed?
  norm_ = ::sqrt(norm_);
  return norm_;
}


double Tensor::ddot(Ref<Tensor>& other) const {
  const blasint unit = 1;
  #ifndef DISK_BASED_SMITH
  sync();
  const double* buffer1 = (double *)        file_->localdata();
  const double* buffer2 = (double *) other->file()->localdata();
  const blasint dsize = file_->localsize() / sizeof(double);
  MPQC_ASSERT(file()->localsize() == other->file()->localsize());
  double ddotproduct = F77_DDOT(&dsize, buffer1, &unit, buffer2, &unit);

  sync();
  Ref<MessageGrp> msg_ = MessageGrp::get_default_messagegrp();
  msg_->sum(ddotproduct);
#else
  file_->clear();
  other->file_->clear();
  double* aux_array = new double[cachesize];
  double* aux_array2 = new double[cachesize];
  double ddotproduct = 0.0;
  long size_back = filesize_; 
  size_t size_now = 0LU;
  while (size_back > 0L) {
    const blasint rsize = min(cachesize, size_back);
    const size_t readsize = rsize * sizeof(double);
    const size_t position = size_now * sizeof(double);
    file_->seekg(position);
    file_->read((char*)aux_array, readsize);
    other->file_->seekg(position);
    other->file_->read((char*)aux_array2, readsize);
    ddotproduct += F77_DDOT(&rsize, aux_array, &unit, aux_array2, &unit);
    size_back -= cachesize;
    size_now += cachesize;
  }
  delete[] aux_array;
  delete[] aux_array2;
#endif
  return ddotproduct;
}

double sc::RMS(const Tensor& t) {
  return t.norm() / t.get_filesize();
}

void Tensor::print(const std::string& label,
                   std::ostream& os) const {
#ifndef DISK_BASED_SMITH 
  os << indent << label << " Tensor" << std::endl;
  typedef std::map<long, long>::const_iterator iter_t;
  for(iter_t i = hash_table_.begin();
      i != hash_table_.end();
      ++i) {

    iter_t ii = i;  ++ii;
    if (ii != hash_table_.end()) {
      long doffset = i->second;
      long dsize   = ii->second-doffset;
      distsize_t offset = (distsize_t)doffset * sizeof(double);
      const int    size =        (int)dsize   * sizeof(double); // in byte
      double* buffer = (double*) file()->obtain_readonly(offset, size);
      os << indent << "tile " << i->first << std::endl;
      for(int k=0; k<dsize; ++k)
        os << indent << "tile[" << k << "] = " << buffer[k] << std::endl;
      file()->release_readonly(buffer, offset, size);
    }

  }
#endif
}

