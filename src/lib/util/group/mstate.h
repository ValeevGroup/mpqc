//
// mstate.h
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
#pragma interface
#endif

#ifndef _util_group_mstate_h
#define _util_group_mstate_h

#include <util/state/state.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/message.h>

namespace sc {

/** The MsgStateSend is an abstract base class that sends objects
    to nodes in a MessageGrp.
*/
class MsgStateSend: public StateOut {
  private:
    // do not allow copy constructor or assignment
    MsgStateSend(const MsgStateSend&);
    void operator=(const MsgStateSend&);
  protected:
    Ref<MessageGrp> grp;
    int nbuf; // the number of bytes used in the buffer
    int bufsize; // the allocated size of the data buffer
    char* buffer; // the data buffer
    char* send_buffer; // the buffer used to send data (includes nbuf)
    int nheader; // nbuf + nheader = the number of bytes in send_buffer to send
    int* nbuf_buffer; // the pointer to the nbuf stored in the buffer

    int put_array_void(const void*, int);
  public:
    MsgStateSend(const Ref<MessageGrp>&);
    virtual ~MsgStateSend();

    /// Specializations must implement flush().
    virtual void flush() = 0;

    /** The buffer size of statein and stateout objects that communicate
        with each other must match. */
    void set_buffer_size(int);

    /** I only need to override put(const ClassDesc*) but C++ will
        hide all of the other put's so I must override everything. */
    int put(const ClassDesc*);
    int put(char r);
    int put(unsigned int r);
    int put(int r);
    int put(unsigned long r);
    int put(long r);
    int put(float r);
    int put(double r);
    int put(const char*,int);
    int put(const int*,int);
    int put(const unsigned int*,int);
    int put(const long*,int);
    int put(const unsigned long*,int);
    int put(const float*,int);
    int put(const double*,int);
};

/** The MsgStateBufRecv is an abstract base class that
    buffers objects sent through a MessageGrp.
*/
class MsgStateBufRecv: public StateIn {
  private:
    // do not allow copy constructor or assignment
    MsgStateBufRecv(const MsgStateBufRecv&);
    void operator=(const MsgStateBufRecv&);
  protected:
    Ref<MessageGrp> grp;
    int nbuf; // the number of bytes used in the buffer
    int ibuf; // the current pointer withing the buffer
    int bufsize; // the allocated size of the buffer
    char* buffer; // the data buffer
    char* send_buffer; // the buffer used to send data (includes nbuf)
    int nheader; // nbuf + nheader = the number of bytes in send_buffer to send
    int* nbuf_buffer; // the pointer to the nbuf stored in the buffer

    int get_array_void(void*,int);

    /// Specializations must implement next_buffer().
    virtual void next_buffer() = 0;
  public:
    /// MsgStateBufRecv can be initialized with a MessageGrp.
    MsgStateBufRecv(const Ref<MessageGrp>&);
    /// Use the default MessageGrp.
    MsgStateBufRecv();

    virtual ~MsgStateBufRecv();

    /** The buffer size of statein and stateout objects that communicate
        with each other must match. */
    void set_buffer_size(int);
};

/** The MsgStateRecv is an abstract base class that receives
    objects from nodes in a MessageGrp. */
class MsgStateRecv: public MsgStateBufRecv {
  private:
    // do not allow copy constructor or assignment
    MsgStateRecv(const MsgStateRecv&);
    void operator=(const MsgStateRecv&);
  public:
    /// MsgStateRecv must be initialized with a MessageGrp.
    MsgStateRecv(const Ref<MessageGrp>&);

    virtual ~MsgStateRecv();

    /** Returns the version of the ClassDesc.  This assumes that
        the version of the remote class is the same as that of
        the local class. */
    int version(const ClassDesc*);

    /** I only need to override get(ClassDesc**) but C++ will hide
        all of the other get's so I must override everything. */
    int get(const ClassDesc**);
    int get(char&r, const char *key = 0);
    int get(unsigned int&r, const char *key = 0);
    int get(int&r, const char *key = 0);
    int get(unsigned long&r, const char *key = 0);
    int get(long&r, const char *key = 0);
    int get(float&r, const char *key = 0);
    int get(double&r, const char *key = 0);
    int get(char*&);
    int get(unsigned int*&);
    int get(int*&);
    int get(unsigned long*&);
    int get(long*&);
    int get(float*&);
    int get(double*&);
};

/** StateSend is a concrete specialization of
    MsgStateSend that does the send part of point to
    point communication in a MessageGrp. */
class StateSend: public MsgStateSend {
  private:
    // do not allow copy constructor or assignment
    StateSend(const StateSend&);
    void operator=(const StateSend&);
  private:
    int type_;
    int target_;
  public:
    /// Create a StateSend given a MessageGrp.
    StateSend(const Ref<MessageGrp>&);

    ~StateSend();
    /// Specify the target node.
    void target(int);
    /// Return the target.
    int get_target() const { return target_; }
    /// Specify the type.
    void type(int);
    /// Return the type.
    int get_type() const { return type_; }
    /// Flush the buffer.
    void flush();
};

/** StateRecv is a concrete specialization of
    MsgStateRecv that does the receive part of point to
    point communication in a MessageGrp. */
class StateRecv: public MsgStateRecv {
  private:
    // do not allow copy constructor or assignment
    StateRecv(const StateRecv&);
    void operator=(const StateRecv&);
  private:
    int source_;
    int type_;
    int last_source_;
    int last_type_;
  protected:
    void next_buffer();
  public:
    /// Create a StateRecv given a MessageGrp.
    StateRecv(const Ref<MessageGrp>&);
    /// Specify the source node.
    void source(int);
    /// Specify the message type.
    void type(int);
    /// Return the source of the last message received.
    int last_source();
    /// Return the type of the last message received.
    int last_type();
};

/** BcastStateSend does the send part of a broadcast of an object
    to all nodes.  Only one node uses a BcastStateSend and the rest
    must use a BcastStateRecv. */
class BcastStateSend: public MsgStateSend {
  private:
    // do not allow copy constructor or assignment
    BcastStateSend(const BcastStateSend&);
    void operator=(const BcastStateSend&);
  public:
    /// Create the BcastStateSend.
    BcastStateSend(const Ref<MessageGrp>&);

    ~BcastStateSend();
    /// Flush the data remaining in the buffer.
    void flush();
};

/** BcastStateRecv does the receive part of a broadcast of an
    object to all nodes.  Only one node uses a BcastStateSend and
    the rest must use a BcastStateRecv. */
class BcastStateRecv: public MsgStateRecv {
  private:
    // do not allow copy constructor or assignment
    BcastStateRecv(const BcastStateRecv&);
    void operator=(const BcastStateRecv&);
  protected:
    int source_;
    void next_buffer();
  public:
    /// Create the BcastStateRecv.
    BcastStateRecv(const Ref<MessageGrp>&, int source = 0);
    /// Set the source node.
    void source(int s);
};

/** This creates and forwards/retrieves data from either a BcastStateRecv
    or a BcastStateSend depending on the value of the argument to
    constructor. */
class BcastState {
  private:
    BcastStateRecv *recv_;
    BcastStateSend *send_;
  public:
    /// Create a BcastState object.  The default source is node 0.
    BcastState(const Ref<MessageGrp> &, int source = 0);

    ~BcastState();

    /** @name Broadcast Members
        Broadcast data to all nodes.  After these are called
        for a group of data the flush member must be called
        to force the source node to actually write the data. */
    //@{
    void bcast(int &);
    void bcast(double &);
    void bcast(int *&, int);
    void bcast(double *&, int);
    template <class T> void bcast(Ref<T>&a)
        {
          if (recv_) {
              a << SavableState::restore_state(*recv_);
            }
          else if (send_) {
              SavableState::save_state(a.pointer(),*send_);
            }
        }
    //@}

    /** Force data to be written.  Data is not otherwise written
        until the buffer is full. */
    void flush();

    /** Call the StateOut or StateIn
        forget_references member. */
    void forget_references();

    /// Controls the amount of data that is buffered before it is sent.
    void set_buffer_size(int);
};

/** BcastStateBin reads a file in written by
    StateInBin on node 0 and broadcasts it to all nodes
    so state can be simultaneously restored on all nodes. */
class BcastStateInBin: public MsgStateBufRecv {
  private:
    // do not allow copy constructor or assignment
    BcastStateInBin(const BcastStateRecv&);
    void operator=(const BcastStateRecv&);
  protected:
    int opened_;
    int file_position_;
    std::streambuf *buf_;

    void next_buffer();
    int get_array_void(void*, int);
  public:
    /// Create the BcastStateRecv using the default MessageGrp.
    BcastStateInBin(const Ref<KeyVal> &);
    /// Create the BcastStateRecv.
    BcastStateInBin(const Ref<MessageGrp>&, const char *filename);

    ~BcastStateInBin();

    virtual int open(const char *name);
    virtual void close();

    void seek(int loc);
    int seekable();
    int tell();
    int use_directory();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
