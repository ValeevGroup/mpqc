//
// actmsg.h
//
// based on: memamsg.h
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

#ifndef _util_group_actmsg_h
#define _util_group_actmsg_h

#include <iostream>

#include <util/group/thread.h>
#include <util/group/memmsg.h>
#include <util/group/mstate.h>

namespace sc {

class ActiveMessageGrp;

/** Derivatives of ActiveMessage can be constructed in one process
    and executed in another by using ActiveMessageGrp.
*/
class ActiveMessage: virtual public SavableState {
  public:
    ActiveMessage() {}
    ActiveMessage(StateIn &s): SavableState(s) {}
    void save_data_state(StateOut &) {}
    /** This is called when ActiveMessageGrp is used to send an
        ActiveMessage object to a process. */
    virtual void run(int sender, int type, ActiveMessageGrp *context) = 0;
};

/// This is an ActiveMessage derivative used for testing. It writes an
/// integer to the output.
class ActiveMessageEcho: public ActiveMessage {
    int i_;
  public:
    ActiveMessageEcho(StateIn &);
    ActiveMessageEcho(int i): i_(i) {}
    void save_data_state(StateOut &);
    void run(int sender, int type, ActiveMessageGrp *context);
};

/// This is a help class that is used by ActiveMessageGrp. It is used to
/// receive and execute ActiveMessage objects.
class ActiveMessageThread: public Thread {
  private:
    ActiveMessageGrp *context_;
    Ref<StateRecv> in_;
    unsigned int nreq_recd_;
  public:
    ActiveMessageThread(const Ref<StateRecv> &,
                        ActiveMessageGrp *context);
    void run();
    int run_one();
    unsigned int nreq_recd() { return nreq_recd_; }
    void set_nreq_recd(unsigned int val) { nreq_recd_ = val; }
};

/** ActiveMessageGrp provides an implemention of active messages that
    sends objects derived from ActiveMessage to remote processes
    and causes their run member to be executed there.
*/
class ActiveMessageGrp : public DescribedClass {
  protected:
    int active_;
    unsigned int *nreq_sent_;
    ActiveMessageThread **thread_;

    int statesend_type_;

    Ref<MessageGrp> msg_;
    Ref<ThreadGrp> thr_;

    void init(const Ref<MessageGrp>& msg,
              const Ref<ThreadGrp>& thr);

  public:
    /** Construct an ActiveMessageGrp using a MessageGrp and a ThreadGrp. */
    ActiveMessageGrp(const Ref<MessageGrp>& msg, const Ref<ThreadGrp>& thr);
    /** A KeyVal CTOR for ActiveMessageGrp. */
    ActiveMessageGrp(const Ref<KeyVal>&);
    ~ActiveMessageGrp();

    /** Each thread using the ActiveMessageGrp needs its own StateSend
        object. This cannot be called concurrently by multiple threads. */
    Ref<StateSend> get_statesend();

    /** Send the active message to node.  The give StateOut must not be
        concurrently used by any other thread.  This member can be called
        concurrently by multiple threads. */
    void send(int node,
              const Ref<StateSend> &,
              const Ref<ActiveMessage> &);

    /** Make the object ready to process messages.  This will also
        synchronizes the nodes.  This must be called by only one thread. */
    void activate();
    /** Processes all outstanding messages and disable processing of
        messages. This will also synchronize the processes.  This must be
        called by only one thread.  */
    void deactivate();
    /** Synchronize all of the processes in this group. */
    void sync();

    /** Return the MessageGrp used to implement the ActiveMessageGrp. */
    Ref<MessageGrp> messagegrp() { return msg_; }

    /// Return the number of processes.
    int n() const { return msg_->n(); }
    /// Return my process identifier, starting at zero.
    int me() const { return msg_->me(); }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
