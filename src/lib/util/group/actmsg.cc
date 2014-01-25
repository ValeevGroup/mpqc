//
// actmsg.cc
//
// based on: memamsg.cc
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

#include "actmsg.h"
#include <util/misc/formio.h>
#include <util/group/mstate.h>

using namespace std;
using namespace sc;

///////////////////////////////////////////////////////////////////////
// The ActiveMessage class

static sc::ClassDesc ActiveMessage_cd(
  typeid(ActiveMessage),"ActiveMessage",1,"virtual public SavableState",
  0, 0, 0);

///////////////////////////////////////////////////////////////////////
// The ActiveMessageEcho class

static sc::ClassDesc ActiveMessageEcho_cd(
  typeid(ActiveMessageEcho),"ActiveMessageEcho",1,"public ActiveMessage",
  0, 0, create<ActiveMessageEcho>);

ActiveMessageEcho::ActiveMessageEcho(StateIn &s):
  SavableState(s),
  ActiveMessage(s)
{
  s.get(i_);
}

void
ActiveMessageEcho::save_data_state(StateOut &s)
{
  ActiveMessage::save_data_state(s);
  s.put(i_);
}

void
ActiveMessageEcho::run(int sender, int type, ActiveMessageGrp *context)
{
  int me = context->messagegrp()->me();
  ExEnv::outn() << " on " << me
                << " got " << i_
                << " from " << sender
                << " type " << type
                << std::endl;
}

///////////////////////////////////////////////////////////////////////
// The ActiveMessageGrp class

ActiveMessageThread::ActiveMessageThread(const Ref<StateRecv> &in,
                                         ActiveMessageGrp *context)
{
  nreq_recd_ = 0;
  in_ = in;
  context_ = context;
}

void
ActiveMessageThread::run()
{
  while (run_one());
}

int
ActiveMessageThread::run_one()
{
  in_->source(MessageGrp::AnySender);
  in_->type(MessageGrp::AnyType);

  int run;
  in_->get(run);
  int type = in_->last_type();
  int source = in_->last_source();

  if (!run) return 0;

  in_->type(type);
  in_->source(source);
  Ref<ActiveMessage> amsg;
  amsg << SavableState::restore_state(*in_.pointer());

  if (amsg == 0) {
      std::cout << "ActiveMessageThread::run_one(): "
                << "got a null ActiveMessage object on "
                << context_->me()
                << std::endl;
      abort();
    }

  amsg->run(source,type,context_);
  nreq_recd_++;

  return 1;
}

///////////////////////////////////////////////////////////////////////
// The ActiveMessageGrp class

ActiveMessageGrp::ActiveMessageGrp(const Ref<MessageGrp>& msg,
                                   const Ref<ThreadGrp>& thr)
{
  // Create an independent communication context.
  init(msg,thr);
}

ActiveMessageGrp::ActiveMessageGrp(const Ref<KeyVal>& keyval)
{
  KeyValValueRefDescribedClass defmsg(MessageGrp::get_default_messagegrp());
  Ref<MessageGrp> msg;
  msg << keyval->describedclassvalue("messagegrp",defmsg);

  KeyValValueRefDescribedClass defthr(ThreadGrp::get_default_threadgrp());
  Ref<ThreadGrp> thr;
  thr <<  keyval->describedclassvalue("threadgrp",defthr);

  init(msg,thr);
}

void
ActiveMessageGrp::init(const Ref<MessageGrp>& msg,
                       const Ref<ThreadGrp>& thr)
{
  statesend_type_ = 0;
  active_ = 0;

  msg_ = msg->clone();
  thr_ = thr->clone(2);

  int nthread = thr_->nthread();
  thread_ = new ActiveMessageThread*[nthread-1];
  thr_->add_thread(0,0);
  for (int i=1; i<nthread; i++) {
      Ref<StateRecv> in = new StateRecv(msg_);
      in->source(-1);
      thread_[i-1] = new ActiveMessageThread(in, this);
      thr_->add_thread(i,thread_[i-1]);
    }

  nreq_sent_ = new unsigned int[n()];
  memset(nreq_sent_, 0, sizeof(unsigned int)*n());
}

ActiveMessageGrp::~ActiveMessageGrp()
{
  deactivate();
  for (int i=0; i<thr_->nthread()-1; i++) {
      delete thread_[i];
    }
  delete[] thread_;
  delete[] nreq_sent_;
}

Ref<StateSend>
ActiveMessageGrp::get_statesend()
{
  Ref<StateSend> out = new StateSend(msg_);
  out->type(statesend_type_++);
  return out;
}

void
ActiveMessageGrp::send(int node, const Ref<StateSend> &out,
                       const Ref<ActiveMessage> &amsg)
{
  if (node == msg_->me()) {
      amsg->run(node,out->get_type(),this);
      return;
    }
  int run = 1;
  out->forget_references();
  out->target(node);
  out->put(run);
  SavableState::save_state(amsg.pointer(),*out.pointer());
  out->flush();
  nreq_sent_[node]++;
}

void
ActiveMessageGrp::activate()
{
//   std::cout << "ActiveMessageGrp::activate() called" << std::endl;

  // Only remote requests require the handler.  There are only remote
  // requests if there is more than one node.
  if (n() == 1) return;

  if (thr_->nthread() < 2) {
      ExEnv::outn() << "ActiveMessageGrp didn't get enough threads" << endl;
      abort();
    }

  if (active_) return;
  active_ = 1;

//   std::cout << "ActiveMessageGrp::activate() starting threads" << std::endl;
  thr_->start_threads();

  msg_->sync();
}

void
ActiveMessageGrp::deactivate()
{
//   std::cout << "ActiveMessageGrp::deactivate() called" << std::endl;

  if (!active_) return;

  active_ = 0;

//   std::cout << "ActiveMessageGrp::deactivate() computing nreq_sent_" << std::endl;

  msg_->sum(nreq_sent_, n());

//   std::cout << "ActiveMessageGrp::deactivate() sending shutdown message" << std::endl;

  // send a shutdown message
  int run = 0;
  Ref<StateSend> out = get_statesend();

  out->target(msg_->me());
  out->put(run);
  out->flush();

  thr_->wait_threads();

//   for (int i=0; i<n(); i++) {
//       if (msg_->me()) {
//           std::cout << "nreq_sent_[" << i << "] = " << nreq_sent_[i] << std::endl;
//         }
//     }

  unsigned int nreq_recd = 0;
  for (int i=0; i<thr_->nthread()-1; i++) {
      nreq_recd += thread_[i]->nreq_recd();
      thread_[i]->set_nreq_recd(0);
    }
  int n_outstanding = nreq_sent_[me()] - nreq_recd;
//   std::cout << "on " << msg_->me() << " processing " << n_outstanding
//             << " outstanding messages "
//             << " received by me " << nreq_recd
//             << " sent to me " << nreq_sent_[me()]
//             << std::endl;
  for (int i=0; i<n_outstanding; i++) {
      thread_[0]->run_one();
    }
  memset(nreq_sent_, 0, sizeof(unsigned int)*n());
  msg_->sync();
}

void
ActiveMessageGrp::sync()
{
  if (active_) {
      deactivate();
      activate();
    }
  else {
      msg_->sync();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
