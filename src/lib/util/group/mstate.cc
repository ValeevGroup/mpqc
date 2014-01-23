//
// mstate.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/group/mstate.h>

#include <util/state/translate.h>

using namespace std;
using namespace sc;

#define DEBUG 0

// This sets up a communication buffer.  It is made up of a of
// an integer that gives the number of bytes used in the buffer
// by the data region of size bufsize.
static
void
obtain_buffer(int*& nbuf_buffer, char*& send_buffer, int& nheader,
              char*& buffer, int& bufsize, int size)
{
  if (size == bufsize) return;
  if (send_buffer) delete[] (int*) send_buffer;

  bufsize = size;

  int min_bytes_to_allocate = bufsize + sizeof(int);
  int ints_to_allocate = min_bytes_to_allocate/sizeof(int);
  if (min_bytes_to_allocate%sizeof(int)) ints_to_allocate++;

  nheader = sizeof(int);
  int * isend_buffer = new int[ints_to_allocate];
  send_buffer = (char*) isend_buffer;
  buffer = (char*) & isend_buffer[1];
  nbuf_buffer = isend_buffer;
}

static
void
release_buffer(char* send_buffer)
{
  if (send_buffer) delete[] (int*)send_buffer;
}

///////////////////////////////////////////////////////////////////////////
// MsgStateSend member functions

MsgStateSend::MsgStateSend(const Ref<MessageGrp>&grp_):
  grp(grp_)
{
  nbuf = 0;
  bufsize = 0;
  send_buffer = 0;
  node_to_node_ = 1;
  obtain_buffer(nbuf_buffer,send_buffer,nheader,buffer,bufsize,8192);
}

MsgStateSend::~MsgStateSend()
{
  release_buffer(send_buffer);
}

void
MsgStateSend::set_buffer_size(int size)
{
  flush();
  obtain_buffer(nbuf_buffer,send_buffer,nheader,buffer,bufsize,size);
}

int
MsgStateSend::put_array_void(const void* vd, int n)
{
  const char* d = (const char*) vd;
  int remaining = n;

  while (remaining) {
      if (nbuf == bufsize) flush();
      int ncurrent;
      if (bufsize - nbuf < remaining) {
          ncurrent = bufsize - nbuf;
        }
      else {
          ncurrent = remaining;
        }
      memcpy(&buffer[nbuf],d,ncurrent);
      remaining -= ncurrent;
      nbuf += ncurrent;
      d = &d[ncurrent];
    }
  return n;
}

int
MsgStateSend::put(const ClassDesc*cd)
{
  int index = grp->classdesc_to_index(cd);
  return StateOut::put(index);
}

int
MsgStateSend::put(const std::string& d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(char d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(unsigned int d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(int d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(unsigned long d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(long d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(bool d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(float d)
{
  return StateOut::put(d);
}


int
MsgStateSend::put(double d)
{
  return StateOut::put(d);
}

int
MsgStateSend::put(const char* d, int n)
{
  return StateOut::put(d, n);
}

int
MsgStateSend::put(const unsigned int* d, int n)
{
  return StateOut::put(d, n);
}

int
MsgStateSend::put(const int* d, int n)
{
  return StateOut::put(d, n);
}

int
MsgStateSend::put(const unsigned long* d, int n)
{
  return StateOut::put(d, n);
}

int
MsgStateSend::put(const long* d, int n)
{
  return StateOut::put(d, n);
}

int
MsgStateSend::put(const float* d, int n)
{
  return StateOut::put(d, n);
}

int
MsgStateSend::put(const double* d, int n)
{
  return StateOut::put(d, n);
}

///////////////////////////////////////////////////////////////////////////
// MsgStateBufRecv member functions

static ClassDesc MsgStateBufRecv_cd(
  typeid(MsgStateBufRecv),"MsgStateBufRecv",1,"public StateIn",
  0, 0, 0);

MsgStateBufRecv::MsgStateBufRecv()
{
  grp = MessageGrp::get_default_messagegrp();
  nbuf = 0;
  ibuf = 0;
  send_buffer = 0;
  bufsize = 0;
  obtain_buffer(nbuf_buffer,send_buffer,nheader,buffer,bufsize,8192);
}

MsgStateBufRecv::MsgStateBufRecv(const Ref<MessageGrp>&grp_):
  grp(grp_)
{
  nbuf = 0;
  ibuf = 0;
  send_buffer = 0;
  bufsize = 0;
  obtain_buffer(nbuf_buffer,send_buffer,nheader,buffer,bufsize,8192);
}

MsgStateBufRecv::~MsgStateBufRecv()
{
  if (ibuf && (nbuf != ibuf)) {
      ExEnv::errn() << scprintf("MsgStateBufRecv::~MsgStateBufRecv(): buffer still has"
              " %d bytes of data on %d\n", nbuf - ibuf, grp->me());
    }
  release_buffer(send_buffer);
}

void
MsgStateBufRecv::set_buffer_size(int size)
{
  if (ibuf && (nbuf != ibuf)) {
      ExEnv::errn() << "MsgStateBufRecv::set_buffer_size(): old buffer has data"
           << endl;
    }
  obtain_buffer(nbuf_buffer, send_buffer, nheader, buffer, bufsize, size);
}

int
MsgStateBufRecv::get_array_void(void* vd, int n)
{
  char* d = (char*) vd;

  int remaining = n;

  while (remaining) {
      if (ibuf == nbuf) next_buffer();
      int ncurrent;
      if (nbuf - ibuf < remaining) {
          ncurrent = nbuf - ibuf;
        }
      else {
          ncurrent = remaining;
        }
      memcpy(d,&buffer[ibuf],ncurrent);
      remaining -= ncurrent;
      ibuf += ncurrent;
      d = &d[ncurrent];
    }

  return n;
}

///////////////////////////////////////////////////////////////////////////
// MsgStateRecv member functions

MsgStateRecv::MsgStateRecv(const Ref<MessageGrp>&grp_):
  MsgStateBufRecv(grp_)
{
  node_to_node_ = 1;
}

MsgStateRecv::~MsgStateRecv()
{
}

int
MsgStateRecv::version(const ClassDesc* cd)
{
  if (!cd) return -1;
  return cd->version();
}

int
MsgStateRecv::get(const ClassDesc**cd)
{
  int index;
  int r = StateIn::get(index);
  *cd = grp->index_to_classdesc(index);
  if (!*cd) {
      ExEnv::errn() << "MsgStateRecvt::get(const ClassDesc**cd): "
           << "class not available on this processor:"
           << endl;
      ExEnv::errn() << " index = " << index << endl;
      abort();
    }
  return r;
}

int
MsgStateRecv::get(std::string& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(char& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(int& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(unsigned int& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(long& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(unsigned long& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(bool& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(float& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(double& d, const char *key)
{
  return StateIn::get(d,key);
}

int
MsgStateRecv::get(char*& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(unsigned int*& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(int*& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(unsigned long*& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(long*& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(float*& d)
{
  return StateIn::get(d);
}

int
MsgStateRecv::get(double*& d)
{
  return StateIn::get(d);
}

///////////////////////////////////////////////////////////////////////////
// StateSend member functions

StateSend::StateSend(const Ref<MessageGrp>&grp_):
  MsgStateSend(grp_),
  type_(0),
  target_(0)
{
}

StateSend::~StateSend()
{
  flush();
}

void
StateSend::flush()
{
  if (nbuf == 0) return;
  *nbuf_buffer = nbuf;
  translate_->translator()->to_external(nbuf_buffer,1);
  grp->raw_sendt(target_, type_, send_buffer, nbuf + nheader);
  nbuf = 0;
}

void
StateSend::target(int t)
{
  target_ = t;
  ps_.clear();
}

void
StateSend::type(int t)
{
  type_ = t;
}

///////////////////////////////////////////////////////////////////////////
// StateRecv member functions

StateRecv::StateRecv(const Ref<MessageGrp>&grp_):
  MsgStateRecv(grp_),
  source_(0),
  type_(-1),
  last_source_(-1),
  last_type_(-1)
{
}

void
StateRecv::next_buffer()
{
  MessageGrp::MessageInfo info;
  grp->raw_recvt(source_, type_, send_buffer, bufsize+nheader,&info);
  translate_->translator()->to_native(nbuf_buffer,1);
  nbuf = *nbuf_buffer;
  ibuf = 0;
  last_type_ = info.type();
  last_source_ = info.sender();
}

void
StateRecv::source(int s)
{
  source_ = s;
  ps_.clear();
}

void
StateRecv::type(int t)
{
  type_ = t;
}

int
StateRecv::last_type()
{
  return last_type_;
}

int
StateRecv::last_source()
{
  return last_source_;
}

///////////////////////////////////////////////////////////////////////////
// BcastStateSend member functions

BcastStateSend::BcastStateSend(const Ref<MessageGrp>&grp_):
  MsgStateSend(grp_)
{
}

BcastStateSend::~BcastStateSend()
{
  flush();
}

void
BcastStateSend::flush()
{
  if (nbuf == 0) return;
  *nbuf_buffer = nbuf;
  translate_->translator()->to_external(nbuf_buffer,1);
  grp->raw_bcast(send_buffer, nbuf + nheader, grp->me());
  nbuf = 0;
}

///////////////////////////////////////////////////////////////////////////
// BcastStateRecv member functions

BcastStateRecv::BcastStateRecv(const Ref<MessageGrp>&grp_, int s):
  MsgStateRecv(grp_)
{
  source(s);
}

void
BcastStateRecv::source(int s)
{
  if (s == grp->me()) {
      ExEnv::errn() << scprintf("BcastStateRecv::source(%d): cannot receive my own"
              " broadcast\n", s);
      abort();
    }
  source_ = s;
  ps_.clear();
}

void
BcastStateRecv::next_buffer()
{
  grp->raw_bcast(send_buffer, bufsize+nheader, source_);
  translate_->translator()->to_native(nbuf_buffer,1);
  nbuf = *nbuf_buffer;
  ibuf = 0;
}

///////////////////////////////////////////////////////////////////////////
// BcastState member functions

BcastState::BcastState(const Ref<MessageGrp> &grp, int source)
{
  if (grp->n() == 1) {
      recv_ = 0;
      send_ = 0;
    }
  else if (grp->me() == source) {
      recv_ = 0;
      send_ = new BcastStateSend(grp);
    }
  else {
      recv_ = new BcastStateRecv(grp,source);
      send_ = 0;
    }
}

BcastState::~BcastState()
{
  delete recv_;
  delete send_;
}

void
BcastState::bcast(int &a)
{
  if (recv_) recv_->get(a);
  else if (send_) send_->put(a);
}

void
BcastState::bcast(double &a)
{
  if (recv_) recv_->get(a);
  else if (send_) send_->put(a);
}

void
BcastState::bcast(int *&a, int n)
{
  if (recv_) recv_->get(a);
  else if (send_) send_->put(a,n);
}

void
BcastState::bcast(double *&a, int n)
{
  if (recv_) recv_->get(a);
  else if (send_) send_->put(a,n);
}

void
BcastState::flush()
{
  if (send_) send_->flush();
}

void
BcastState::set_buffer_size(int n)
{
  if (send_) send_->set_buffer_size(n);
  if (recv_) recv_->set_buffer_size(n);
}

void
BcastState::forget_references()
{
  if (send_) send_->forget_references();
}

///////////////////////////////////////////////////////////////////////////
// BcastStateRecv member functions

static ClassDesc BcastStateInBin_cd(
  typeid(BcastStateInBin),"BcastStateInBin",1,"public MsgStateBufRecv",
  0, create<BcastStateInBin>, 0);

BcastStateInBin::BcastStateInBin(const Ref<MessageGrp>&grp_,
                                 const char *filename):
  MsgStateBufRecv(grp_)
{
  opened_ = 0;
  open(filename);
}

BcastStateInBin::BcastStateInBin(const Ref<KeyVal> &keyval)
{
  std::string path = keyval->stringvalue("file");
  if (path.empty()) {
      ExEnv::errn() << "StateInBin(const Ref<KeyVal>&): no path given" << endl;
    }
  opened_ = 0;
  open(path.c_str());
}

BcastStateInBin::~BcastStateInBin()
{
  close();
}

void
BcastStateInBin::next_buffer()
{
  if (grp->me() == 0) {
      // fill the buffer
      *nbuf_buffer = buf_->sgetn(buffer,bufsize);
      if (*nbuf_buffer == 0) {
          ExEnv::errn() << "BcastStateInBin: read failed" << endl;
          abort();
        }
      translate_->translator()->to_external(nbuf_buffer,1);
    }
  grp->raw_bcast(send_buffer, bufsize+nheader);
  translate_->translator()->to_native(nbuf_buffer,1);
  nbuf = *nbuf_buffer;
  ibuf = 0;
}

void
BcastStateInBin::close()
{
  if(opened_) delete buf_;
  opened_=0; buf_=0;
  nbuf = 0;
  ibuf = 0;

  classidmap_.clear();
  nextclassid_ = 0;
  classdatamap_.clear();
  ps_.clear();
}

int
BcastStateInBin::open(const char *path)
{
  file_position_ = 0;

  if (grp->me() == 0) {
      if (opened_) close();

      filebuf *fbuf = new filebuf();
      fbuf->open(path, ios::in);
      if (!fbuf->is_open()) {
          ExEnv::errn() << "ERROR: BcastStateInBin: problems opening " << path << endl;
          abort();
        }
      buf_ = fbuf;
      opened_ = 1;
    }

  nbuf = 0;
  ibuf = 0;

  get_header();
  find_and_get_directory();

  return 0;
}

int
BcastStateInBin::tell()
{
  return file_position_;
}

void
BcastStateInBin::seek(int loc)
{
  file_position_ = loc;
#if defined(HAVE_PUBSEEKOFF)
  if (grp->me() == 0) {
      buf_->pubseekoff(loc,ios::beg,ios::in);
#  if  DEBUG
      ExEnv::outn() << "pubseekoff to " << loc << endl;
#  endif
    }
#elif defined(HAVE_SEEKOFF)
  if (grp->me() == 0) {
      buf_->seekoff(loc,ios::beg,ios::in);
#  if  DEBUG
      ExEnv::outn() << "seekoff to " << loc << endl;
#  endif
    }
#endif
  nbuf = 0;
  ibuf = 0;
}

int
BcastStateInBin::seekable()
{
#if defined(HAVE_PUBSEEKOFF) || defined(HAVE_SEEKOFF)
  return 1;
#else
  return 0;
#endif
}

int
BcastStateInBin::use_directory()
{
  return seekable();
}

int
BcastStateInBin::get_array_void(void* vd, int n)
{
  MsgStateBufRecv::get_array_void(vd, n);
  file_position_ += n;
#if DEBUG
  ExEnv::outn() << "Read " << n << " bytes:";
  for (int i=0; i<n; i++) {
      ExEnv::outn() << " " << (int) ((unsigned char*)vd)[i];
    }
  ExEnv::outn() << endl;
#endif
  return n;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
