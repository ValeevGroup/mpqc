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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/group/mstate.h>

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

MsgStateSend::MsgStateSend(const RefMessageGrp&grp_):
  grp(grp_)
{
  nbuf = 0;
  bufsize = 0;
  send_buffer = 0;
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
  return StateOutBinXDR::put(grp->classdesc_to_index(cd));
}

int
MsgStateSend::put(char d)
{
  return StateOutBinXDR::put(d);
}

int
MsgStateSend::put(int d)
{
  return StateOutBinXDR::put(d);
}

int
MsgStateSend::put(float d)
{
  return StateOutBinXDR::put(d);
}


int
MsgStateSend::put(double d)
{
  return StateOutBinXDR::put(d);
}

int
MsgStateSend::put(char* d, int n)
{
  return StateOutBinXDR::put(d, n);
}

int
MsgStateSend::put(int* d, int n)
{
  return StateOutBinXDR::put(d, n);
}

int
MsgStateSend::put(float* d, int n)
{
  return StateOutBinXDR::put(d, n);
}

int
MsgStateSend::put(double* d, int n)
{
  return StateOutBinXDR::put(d, n);
}

///////////////////////////////////////////////////////////////////////////
// MsgStateRecv member functions

MsgStateRecv::MsgStateRecv(const RefMessageGrp&grp_):
  grp(grp_)
{
  nbuf = 0;
  ibuf = 0;
  send_buffer = 0;
  bufsize = 0;
  obtain_buffer(nbuf_buffer,send_buffer,nheader,buffer,bufsize,8192);
}

MsgStateRecv::~MsgStateRecv()
{
  if (ibuf && (nbuf != ibuf)) {
      cerr << scprintf("MsgStateRecv::~MsgStateRecv(): buffer still has"
              " %d bytes of data\n", nbuf - ibuf);
    }
  release_buffer(send_buffer);
}

void
MsgStateRecv::set_buffer_size(int size)
{
  if (ibuf && (nbuf != ibuf)) {
      cerr << scprintf("MsgStateRecv::set_buffer_size(): old buffer has data\n");
    }
  obtain_buffer(nbuf_buffer, send_buffer, nheader, buffer, bufsize, size);
}

int
MsgStateRecv::get_array_void(void* vd, int n)
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

int
MsgStateRecv::get(const ClassDesc**cd)
{
  int index;
  int r = StateInBinXDR::get(index);
  *cd = grp->index_to_classdesc(index);
  if (!*cd) {
      cerr << scprintf("MsgStateRecvt::get(const ClassDesc**cd): "
              "class not available on this processor\n");
      abort();
    }
  return r;
}

int
MsgStateRecv::get(char& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(int& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(float& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(double& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(char*& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(int*& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(float*& d)
{
  return StateInBinXDR::get(d);
}

int
MsgStateRecv::get(double*& d)
{
  return StateInBinXDR::get(d);
}

///////////////////////////////////////////////////////////////////////////
// StateSend member functions

StateSend::StateSend(const RefMessageGrp&grp_):
  MsgStateSend(grp_),
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
  translate(nbuf_buffer);
  grp->raw_send(target_, send_buffer, nbuf + nheader);
  nbuf = 0;
}

void
StateSend::target(int t)
{
  target_ = t;
  forget_references();
}

///////////////////////////////////////////////////////////////////////////
// StateRecv member functions

StateRecv::StateRecv(const RefMessageGrp&grp_):
  MsgStateRecv(grp_),
  source_(0)
{
}

void
StateRecv::next_buffer()
{
  grp->raw_recv(source_, send_buffer, bufsize+nheader);
  translate(nbuf_buffer);
  nbuf = *nbuf_buffer;
  ibuf = 0;
}

void
StateRecv::source(int s)
{
  source_ = s;
  forget_references();
}

///////////////////////////////////////////////////////////////////////////
// BcastStateSend member functions

BcastStateSend::BcastStateSend(const RefMessageGrp&grp_):
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
  translate(nbuf_buffer);
  grp->raw_bcast(send_buffer, nbuf + nheader, grp->me());
  nbuf = 0;
}

///////////////////////////////////////////////////////////////////////////
// BcastStateRecv member functions

BcastStateRecv::BcastStateRecv(const RefMessageGrp&grp_, int s):
  MsgStateRecv(grp_)
{
  source(s);
}

void
BcastStateRecv::source(int s)
{
  if (s == grp->me()) {
      cerr << scprintf("BcastStateRecv::source(%d): cannot receive my own"
              " broadcast\n", s);
      abort();
    }
  source_ = s;
  forget_references();
}

void
BcastStateRecv::next_buffer()
{
  grp->raw_bcast(send_buffer, bufsize+nheader, source_);
  translate(nbuf_buffer);
  nbuf = *nbuf_buffer;
  ibuf = 0;
}

///////////////////////////////////////////////////////////////////////////
// BcastState member functions

BcastState::BcastState(const RefMessageGrp &grp, int source)
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
BcastState::bcast(SSRefBase &a)
{
  if (recv_) a.restore_state(*recv_);
  else if (send_) {
      a.save_state(*send_);
    }
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
  if (recv_) recv_->forget_references();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
