//
// message.h
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

#ifndef _util_group_message_h
#define _util_group_message_h

#include <map>
#include <set>

#include <math.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <util/group/topology.h>

namespace sc {

template <class T>
class GrpReduce {
  public:
    GrpReduce() {}
    virtual ~GrpReduce() {};
    virtual void reduce(T*target, T*data, int n) = 0;
};

template <class T>
class GrpSumReduce: public GrpReduce<T> {
  public:
    GrpSumReduce() : GrpReduce<T>() {}
    ~GrpSumReduce() {};
    void reduce(T*target, T*data, int nelement);
};

template <class T>
class GrpMinReduce: public GrpReduce<T> {
  public:
    GrpMinReduce() : GrpReduce<T>() {}
    ~GrpMinReduce() {};
    void reduce(T*target, T*data, int nelement);
};

template <class T>
class GrpMaxReduce: public GrpReduce<T> {
  public:
    GrpMaxReduce() : GrpReduce<T>() {}
    ~GrpMaxReduce() {};
    void reduce(T*target, T*data, int nelement);
};

template <class T, class BinaryPredicate = std::less<T> >
class GrpCompareReduce: public GrpReduce<T> {
  public:
    GrpCompareReduce() : GrpReduce<T>() {}
    ~GrpCompareReduce() {};
    void reduce(T*target, T*data, int nelement);
  private:
    BinaryPredicate Op;
};

template <class T>
class GrpArithmeticAndReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement);
};

template <class T>
class GrpArithmeticOrReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement);
};

template <class T>
class GrpArithmeticXOrReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement);
};

template <class T>
class GrpProductReduce: public GrpReduce<T> {
  public:
    void reduce(T*target, T*data, int nelement);
};

template <class T>
class GrpFunctionReduce: public GrpReduce<T> {
  private:
    void (*func_)(T*target,T*data,int nelement);
  public:
    GrpFunctionReduce(void(*func)(T*,T*,int)):func_(func) {}
    void reduce(T*target, T*data, int nelement);
};

/** The MessageGrp abstract class provides
    a mechanism for moving data and objects between
    nodes in a parallel machine. */
class MessageGrp: public DescribedClass {
  public:
    class MessageInfo {
        friend class MessageGrp;
      private:
        int sender_;
        int type_;
        int nbyte_;
      public:
        int sender() const { return sender_; }
        int type() const { return type_; }
        int nbyte() const { return nbyte_; }
    };
    enum { AnyType = -1 };
    enum { AnySender = -1 };
    class MessageHandle {
        friend class MessageGrp;
      private:
        void *id_;
      public:
        MessageHandle(): id_(0) {}
        MessageHandle(const MessageHandle &h): id_(h.id_) {}
    };
  private:
    // These are initialized by the initialize() member (see below).
    int me_;
    int n_;
    int nclass_;
    int gop_max_;
    std::map<ClassDescP,int> classdesc_to_index_;
    ClassDescP *index_to_classdesc_;
  protected:
    /** The classdesc_to_index_ and index_to_classdesc_ arrays
        cannot be initialized by the MessageGrp CTOR, because
        the MessageGrp specialization has not yet been initialized
        and communication is not available.  CTOR's of specializations
        of MessageGrp must call the initialize member in their body
        to complete the initialization process. */
    void initialize(int me, int n);

    Ref<MachineTopology> topology_;

    int debug_;

    void set_sender(MessageInfo *info,int sender) {
      if (info) info->sender_ = sender;
    }
    void set_type(MessageInfo *info,int type) {
      if (info) info->type_ = type;
    }
    void set_nbyte(MessageInfo *info,int nbyte) {
      if (info) info->nbyte_ = nbyte;
    }

    void set_id(MessageHandle *handle,void *id) {
      handle->id_ = id;
    }
    void *get_id(const MessageHandle *handle) {
      return handle->id_;
    }
  public:

    MessageGrp();
    MessageGrp(const Ref<KeyVal>&);

    /** Destroy this MessageGrp. This must be called by all nodes
        participating in this MessageGrp, because it can be a
        collective operation. Note, however, that it cannot be assumed
        to always be collective. */
    virtual ~MessageGrp();
    
    /// Returns the number of processors.
    int n() { return n_; }
    /// Returns my processor number.  In the range [0,n()).
    int me() { return me_; }

    /** Returns a copy of this MessageGrp specialization that provides
        an independent communication context. */
    virtual Ref<MessageGrp> clone(void)=0;

    /** Returns MessageGrp objects that are a subset of this
        MessageGrp. This call is collective--all nodes in this
        MessageGrp must call it. If the commkey is less than zero,
        then a nil pointer will be returned. Otherwise, all nodes
        calling with the same commkey will be put into the same return
        MessageGrp. The rankkey argument is used to order the nodes in
        the new MessageGrp. Nodes with the same rankkey are ordered by
        rank. With the default arguments, split behaves the same as
        clone. */
    virtual Ref<MessageGrp> split(int grpkey=0, int rankkey=0) = 0;

    /** Returns MessageGrp objects that are a subset of this
        MessageGrp. This call is collective--all nodes in this
        MessageGrp must call it. The ordering of the ranks matches
        the ordering of the original ranks.
    */
    virtual Ref<MessageGrp> subset(const std::set<int> &) = 0;
    
    /** The default message group contains the primary message group to
        be used by an application. */
    static void set_default_messagegrp(const Ref<MessageGrp>&);
    /// Returns the default message group.
    static MessageGrp* get_default_messagegrp();

    /** Create a message group.  This routine looks for a -messagegrp
        argument, then the environmental variable MESSAGEGRP to decide which
        specialization of MessageGrp would be appropriate.  The
        argument to -messagegrp should be either string for a
        ParsedKeyVal constructor or a classname.  If this returns
        null, it is up to the programmer to create a MessageGrp. */
    static MessageGrp* initial_messagegrp(int &argc, char** &argv);

    /** @name Send Members

        @param target the recipent node number.
        @param data the data.
        @param ndata the number of data.
    */
    //@{
    virtual void send(int target, const double* data, int ndata);
    virtual void send(int target, const unsigned int* data, int ndata);
    virtual void send(int target, const int* data, int ndata);
    virtual void send(int target, const char* data, int nbyte);
    virtual void send(int target, const unsigned char* data, int nbyte);
    virtual void send(int target, const signed char* data, int nbyte);
    virtual void send(int target, const short* data, int ndata);
    virtual void send(int target, const long* data, int ndata);
    virtual void send(int target, const float* data, int ndata);
    /// This sends a single double datum.
    void send(int target, double data) { send(target,&data,1); }
    /// This sends a single integer datum.
    void send(int target, int data) { send(target,&data,1); }
    virtual void raw_send(int target, const void* data, int nbyte) = 0;
    //@}

    /** @name Send Members for Typed Messages

        @param target the recipent node number.
        @param type the user-chosen type which is used to match the receive.
        @param data the data.
        @param ndata the number of data.
    */
    //@{
    virtual void sendt(int target, int type, const double* data, int ndata,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const unsigned int* data, int ndata,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const int* data, int ndata,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const char* data, int nbyte,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const unsigned char* data, int nbyte,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const signed char* data, int nbyte,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const short* data, int ndata,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const long* data, int ndata,
                       bool rcvrdy=false);
    virtual void sendt(int target, int type, const float* data, int ndata,
                       bool rcvrdy=false);
    /// This sends a single double datum.
    void sendt(int target, int type, double data,
               bool rcvrdy=false) {sendt(target,type,&data,1,rcvrdy);}
    /// This sends a single integer datum.
    void sendt(int target, int type, int data,
               bool rcvrdy=false) {sendt(target,type,&data,1,rcvrdy);}
    virtual void raw_sendt(int target, int type, const void* data, int nbyte,
                           bool rcvrdy=false) = 0;
    //@}

    /** @name Receive Members
        @param sender the sending node number (if -1, any sender will match).
        @param data the data.
        @param ndata the number of data.
     */
    //@{
    virtual void recv(int sender, double* data, int ndata);
    virtual void recv(int sender, unsigned int* data, int ndata);
    virtual void recv(int sender, int* data, int ndata);
    virtual void recv(int sender, char* data, int nbyte);
    virtual void recv(int sender, unsigned char* data, int nbyte);
    virtual void recv(int sender, signed char* data, int nbyte);
    virtual void recv(int sender, short* data, int ndata);
    virtual void recv(int sender, long* data, int ndata);
    virtual void recv(int sender, float* data, int ndata);
    /// This receives a single double datum.
    void recv(int sender, double& data) { recv(sender,&data,1); }
    /// This receives a single integer datum.
    void recv(int sender, int& data) { recv(sender,&data,1); }
    virtual void raw_recv(int sender, void* data, int nbyte,
                          MessageInfo *info=0) = 0;
    //@}

    /** @name Receive Members for Typed Messages
        @param sender the sending node number (if -1, any sender will match).
        @param type the user-chosen type which is used to match the send.
        @param data the data.
        @param ndata the number of data.
     */
    //@{
    virtual void recvt(int sender, int type, double* data, int ndata);
    virtual void recvt(int sender, int type, unsigned int* data, int ndata);
    virtual void recvt(int sender, int type, int* data, int ndata);
    virtual void recvt(int sender, int type, char* data, int nbyte);
    virtual void recvt(int sender, int type, unsigned char* data, int nbyte);
    virtual void recvt(int sender, int type, signed char* data, int nbyte);
    virtual void recvt(int sender, int type, short* data, int ndata);
    virtual void recvt(int sender, int type, long* data, int ndata);
    virtual void recvt(int sender, int type, float* data, int ndata);
    /// This receives a single double datum.
    void recvt(int sender, int type, double& data) {
      recvt(sender,type,&data,1);
    }
    /// This receives a single integer datum.
    void recvt(int sender, int type, int& data) {
      recvt(sender,type,&data,1);
    }
    virtual void raw_recvt(int sender, int type, void* data, int nbyte,
                           MessageInfo *info=0) = 0;
    //@}

    /** @name Non-blocking Send Members
        These send routines are nonblocking. The operation
        is not complete and the \p data buffer cannot be reused
        until a wait completes on the \p handle.

        @param target the recipent node number.
        @param type the user-chosen type which is used to match the receive.
        @param data the data.
        @param ndata the number of data.
        @param handle set to the handle which must be used with wait
        to check for completion.
        @param rcvrdy if true, the receive is guaranteed to have been issued.
        The default is false.

        \sa wait
     */
    //@{
    virtual void nb_sendt(int target, int type,
                          const double* data, int ndata,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const unsigned int* data, int ndata,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const int* data, int ndata,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const char* data, int nbyte,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const unsigned char* data, int nbyte,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const signed char* data, int nbyte,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const short* data, int ndata,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const long* data, int ndata,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    virtual void nb_sendt(int target, int type,
                          const float* data, int ndata,
                          MessageHandle&handle,
                          bool rcvrdy=false);
    void nb_sendt(int target, int type, double data,
                  MessageHandle&handle,
                  bool rcvrdy=false) {
      nb_sendt(target,type,&data,1,handle,rcvrdy);
    }
    void nb_sendt(int target, int type, int data,
                  MessageHandle&handle,
                  bool rcvrdy=false) {
      nb_sendt(target,type,&data,1,handle,rcvrdy);
    }
    virtual void raw_nb_sendt(int target, int type,
                              const void* data, int nbyte,
                              MessageHandle&,
                              bool rcvrdy=false) = 0;
    //@}

    /** @name Non-blocking Receive Members
        These receive routines are nonblocking. The operation
        is not complete and the \p data buffer cannot be used
        until a wait completes on the \p handle.

        @param sender the sending node number (if -1, any sender will match).
        @param type the user-chosen type which is used to match the send.
        @param data the data.
        @param ndata the number of data.
        @param handle set to the handle which must be used with wait
        to check for completion.

        \sa wait
     */
    //@{
    virtual void nb_recvt(int sender, int type, double* data, int ndata,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, unsigned int* data, int ndata,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, int* data, int ndata,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, char* data, int nbyte,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, unsigned char* data, int nbyte,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, signed char* data, int nbyte,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, short* data, int ndata,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, long* data, int ndata,
                          MessageHandle&handle);
    virtual void nb_recvt(int sender, int type, float* data, int ndata,
                          MessageHandle&handle);
    /// This receives a single double datum.
    void nb_recvt(int sender, int type, double& data,
                  MessageHandle&handle) {
      nb_recvt(sender,type,&data,1,handle);
    }
    /// This receives a single integer datum.
    void nb_recvt(int sender, int type, int& data,
                  MessageHandle&handle) {
      nb_recvt(sender,type,&data,1,handle);
    }
    virtual void raw_nb_recvt(int sender, int type,
                              void* data, int nbyte,
                              MessageHandle&) = 0;
    //@}

    /** Wait for an operation to complete.
        @param handle the handle initialized by, for example, nb_recvt.
        @param info if non-null, this will be filled in with information
        about the operation that completed.
    */
    virtual void wait(const MessageHandle&handle,
                      MessageInfo *info=0) = 0;

    /// Ask if a given typed message has been received.
    virtual int probet(int sender, int type, MessageInfo*info=0) = 0;

    /** @name Broadcast Members
        Do broadcasts of various types of data. */
    //@{
    virtual void bcast(double* data, int ndata, int from = 0);
    virtual void bcast(unsigned int* data, int ndata, int from = 0);
    virtual void bcast(int* data, int ndata, int from = 0);
    virtual void bcast(char* data, int nbyte, int from = 0);
    virtual void bcast(unsigned char* data, int nbyte, int from = 0);
    virtual void bcast(signed char* data, int nbyte, int from = 0);
    virtual void bcast(short* data, int ndata, int from = 0);
    virtual void bcast(long* data, int ndata, int from = 0);
    virtual void bcast(float* data, int ndata, int from = 0);
    virtual void raw_bcast(void* data, int nbyte, int from = 0);
    void bcast(double& data, int from = 0) { bcast(&data, 1, from); }
    void bcast(int& data, int from = 0) { bcast(&data, 1, from); }
    //@}

    /** @name Data Collection Members
        Collect data distributed on the nodes to a big array replicated
        on each node. */
    //@{
    virtual void raw_collect(const void *part, const int *lengths,
                             void *whole, int bytes_per_datum=1);
    void collect(const double *part, const int *lengths, double *whole);
    //@}

    /** @name Global Sum Reduction Members */
    //@{
    virtual void sum(double* data, int n, double* = 0, int target = -1);
    virtual void sum(unsigned int* data, int n, unsigned int* = 0, int target = -1);
    virtual void sum(int* data, int n, int* = 0, int target = -1);
    virtual void sum(long* data, int n, long* = 0, int target = -1);
    virtual void sum(unsigned long* data, int n, unsigned long* = 0, int target = -1);
    virtual void sum(char* data, int n, char* = 0, int target = -1);
    virtual void sum(unsigned char* data, int n,
                     unsigned char* = 0, int target = -1);
    virtual void sum(signed char* data, int n,
                     signed char* = 0, int target = -1);
    void sum(double& data) { sum(&data, 1); }
    void sum(int& data) { sum(&data, 1); }
    //@}

    /** @name Global Maximization Members */
    //@{
    virtual void max(double* data, int n, double* = 0, int target = -1);
    virtual void max(int* data, int n, int* = 0, int target = -1);
    virtual void max(unsigned int* data, int n, unsigned int* = 0, int target = -1);
    virtual void max(char* data, int n, char* = 0, int target = -1);
    virtual void max(unsigned char* data, int n,
                     unsigned char* = 0, int target = -1);
    virtual void max(signed char* data, int n,
                     signed char* = 0, int target = -1);
    void max(double& data) { max(&data, 1); }
    void max(int& data) { max(&data, 1); }
    //@}

    /** @name Global Minimization Members */
    //@{
    virtual void min(double* data, int n, double* = 0, int target = -1);
    virtual void min(int* data, int n, int* = 0, int target = -1);
    virtual void min(unsigned int* data, int n, unsigned int* = 0, int target = -1);
    virtual void min(char* data, int n, char* = 0, int target = -1);
    virtual void min(unsigned char* data, int n,
                     unsigned char* = 0, int target = -1);
    virtual void min(signed char* data, int n,
                     signed char* = 0, int target = -1);
    void min(double& data) { min(&data, 1); }
    void min(int& data) { min(&data, 1); }
    //@}

    /** @name Global Generic Reduction Members. */
    //@{
    /// T must be a POD (plain old data) type so that it can be copied bytewise
    template <typename T>
    void reduce(T*, int n, GrpReduce<T>&,
                T*scratch = 0, int target = -1);
    virtual void reduce(double*, int n, GrpReduce<double>&,
                        double*scratch = 0, int target = -1);
    virtual void reduce(int*, int n, GrpReduce<int>&,
                        int*scratch = 0, int target = -1);
    virtual void reduce(unsigned int*, int n, GrpReduce<unsigned int>&,
                        unsigned int*scratch = 0, int target = -1);
    virtual void reduce(char*, int n, GrpReduce<char>&,
                        char*scratch = 0, int target = -1);
    virtual void reduce(unsigned char*, int n, GrpReduce<unsigned char>&,
                        unsigned char*scratch = 0, int target = -1);
    virtual void reduce(signed char*, int n, GrpReduce<signed char>&,
                        signed char*scratch = 0, int target = -1);
    virtual void reduce(short*, int n, GrpReduce<short>&,
                        short*scratch = 0, int target = -1);
    virtual void reduce(float*, int n, GrpReduce<float>&,
                        float*scratch = 0, int target = -1);
    virtual void reduce(long*, int n, GrpReduce<long>&,
                        long*scratch = 0, int target = -1);
    void reduce(double& data, GrpReduce<double>& r) { reduce(&data, 1, r); }
    void reduce(int& data, GrpReduce<int>& r) { reduce(&data, 1, r); }
    //@}

    /// Synchronize all of the processors.
    virtual void sync();

    /// Return the MachineTopology object.
    Ref<MachineTopology> topology() { return topology_; }

    /** @name Global Class Information Members
        Each message group maintains an association of ClassDesc with
        a global index so SavableState information can be sent between
        nodes without needing to send the classname and look up the
        ClassDesc with each transfer.  These routines return information
        about that mapping. */
    //@{
    int classdesc_to_index(const ClassDesc*);
    const ClassDesc* index_to_classdesc(int);
    int nclass() const { return nclass_; }
    //@}
};

struct message_struct {
  void *buf;
  int size;
  int type;
  struct message_struct *p;
  };
typedef struct message_struct message_t;


/** ProcMessageGrp provides a concrete specialization
    of MessageGrp that supports only one node. */
class ProcMessageGrp: public MessageGrp {
  private:
    // Messages are stored in these linked lists
    message_t *sync_messages;
    message_t *type_messages;

    void sendit(message_t *& messages, int dest, int msgtype, const void* buf, int bytes);
    void recvit(message_t *& messages, int source, int type, void* buf, int bytes,
                int& last_size, int& last_type);
        
  public:
    ProcMessageGrp();
    ProcMessageGrp(const Ref<KeyVal>&);
    ~ProcMessageGrp();

    Ref<MessageGrp> clone(void);
    Ref<MessageGrp> split(int grpkey=0, int rankkey=0);
    Ref<MessageGrp> subset(const std::set<int> &);
    
    void raw_send(int target, const void* data, int nbyte);
    void raw_sendt(int target, int type, const void* data, int nbyte,
                   bool rcvrdy=false);
    void raw_recv(int sender, void* data, int nbyte,
                  MessageInfo *info=0);
    void raw_recvt(int sender, int type, void* data, int nbyte,
                   MessageInfo *info=0);
    void raw_bcast(void* data, int nbyte, int from);

    void raw_nb_sendt(int sender, int type,
                      const void* data, int nbyte,
                      MessageHandle&,
                      bool rcvrdy=false);
    void raw_nb_recvt(int sender, int type,
                      void* data, int nbyte,
                      MessageHandle&);
    void wait(const MessageHandle&,
              MessageInfo *info=0);

    int probet(int sender, int type, MessageInfo *info=0);
    void sync();
};

}

#include <util/group/messaget.h>

#endif


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
