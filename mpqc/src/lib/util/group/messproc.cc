//
// messproc.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <util/misc/formio.h>
#include <util/group/message.h>

#define CLASSNAME ProcMessageGrp
#define PARENTS public MessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
ProcMessageGrp::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] =  MessageGrp::_castdown(cd);
  return do_castdowns(casts,cd);
}

ProcMessageGrp::ProcMessageGrp(const RefKeyVal& keyval):
  MessageGrp(keyval)
{
  initialize(0,1);
}

ProcMessageGrp::ProcMessageGrp()
{
  initialize(0,1);
}

ProcMessageGrp::~ProcMessageGrp()
{
}

struct message_struct {
  void *buf;
  int size;
  int type;
  struct message_struct *p;
  };
typedef struct message_struct message_t;

// Messages are stored in these linked lists.
static message_t *sync_messages=0;
static message_t *type_messages=0;

static
void
sendit(message_t *& messages, int dest, int msgtype, void* buf, int bytes)
{
  message_t *msg;
  message_t *I;

  if (dest != 0) {
      cerr << scprintf("messproc.cc:sendit: can only send to 0\n");
      abort();
    }

  msg = (message_t *) malloc(sizeof(message_t));
  if (msg) msg->buf = (char *) malloc(bytes);
  if (!msg || !msg->buf) {
      cerr << scprintf("messproc.cc:sendit: allocation failed\n");
      abort();
    }

  // Put msg at the end of the linked list, because of some bad
  // assumptions made by the mpscf program and libraries.
  msg->p = 0;
  if (!messages) {
      messages = msg;
    }
  else {
      for (I=messages; I->p != 0; I=I->p);
      I->p = msg;
    }

  memcpy(msg->buf,buf,bytes);
  msg->type = msgtype;
  msg->size = bytes;
}

static
void
recvit(message_t *& messages, int source, int type, void* buf, int bytes,
       int& last_size, int& last_type)
{
  message_t *i;
  message_t *last;

  last = 0;
  for (i=messages; i!=0; i = i->p) {
    if (i->type == type || type == -1) {
      if (i->size > bytes) {
        cerr << scprintf(
                "messproc.cc:recvit: message buffer isn't big enough\n");
        abort();
        }
      memcpy(buf,i->buf,i->size);
      last_size = i->size;
      last_type = i->type;

      // Remove the message from the list.
      if (last) {
          last->p = i->p;
        }
      else {
          messages = messages->p;
        }
      free(i->buf);
      free(i);

      return;
      }
    last = i;
    }

  cerr << scprintf(
          "messproc.cc:recvit: tried to receive something that isn't there\n");
  cerr << scprintf("messproc:recvit: tried %d bytes of type %d, ",bytes,type);
  abort();
}

void
ProcMessageGrp::raw_send(int target, void* data, int nbyte)
{
  sendit(sync_messages, target, -1, data, nbyte);
}

void
ProcMessageGrp::raw_sendt(int target, int type, void* data, int nbyte)
{
  sendit(type_messages, target, type, data, nbyte);
}

void
ProcMessageGrp::raw_recv(int sender, void* data, int nbyte)
{
  int last_size, last_type;
  recvit(sync_messages, sender, -1, data, nbyte, last_size, last_type);
  set_last_size(last_size);
  set_last_type(last_type);
  set_last_source(0);
}

void
ProcMessageGrp::raw_recvt(int type, void* data, int nbyte)
{
  int last_size, last_type;
  recvit(type_messages, -1, type, data, nbyte, last_size, last_type);
  set_last_size(last_size);
  set_last_type(last_type);
  set_last_source(0);
}

void
ProcMessageGrp::raw_bcast(void* data, int nbyte, int from)
{
}

int
ProcMessageGrp::probet(int type)
{
  message_t *i;

  for (i=type_messages; i!=0; i = i->p) {
      if (i->type == type || type == -1) {
          set_last_source(0);
          set_last_size(i->size);
          set_last_type(i->type);
          return 1;
        }
    }

  return 0;
}

void
ProcMessageGrp::sync()
{
}

int
ProcMessageGrp::last_source()
{
  return last_source_;
}

int
ProcMessageGrp::last_size()
{
  return last_size_;
}

int
ProcMessageGrp::last_type()
{
  return last_type_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
