
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_mstate_h
#define _util_group_mstate_h

#include <util/state/state.h>
#include <util/group/message.h>

class MsgStateSend: public StateOutXDR {
  private:
    // do not allow copy constructor or assignment
    MsgStateSend(const MsgStateSend&);
    operator=(const MsgStateSend&);
  protected:
    MessageGrp& grp;
    int nbuf; // the number of bytes used in the buffer
    int bufsize; // the allocated size of the data buffer
    char* buffer; // the data buffer
    char* send_buffer; // the buffer used to send data (includes nbuf)
    int nheader; // nbuf + nheader = the number of bytes in send_buffer to send
    int* nbuf_buffer; // the pointer to the nbuf stored in the buffer

    put_array_void(const void*, int);
  public:
    MsgStateSend(MessageGrp&);
    virtual ~MsgStateSend();
    virtual void flush() = 0;

    // The buffer size of statein and stateout objects that
    // communicate with each other must match.
    void set_buffer_size(int);

    // I only need to override put(const ClassDesc*) but
    // C++ will hide all of the other put's so I must
    // override everything.
    int put(const ClassDesc*);
    int put(char r);
    int put(int r);
    int put(float r);
    int put(double r);
    int put(char*,int);
    int put(int*,int);
    int put(float*,int);
    int put(double*,int);
};

class MsgStateRecv: public StateInXDR {
  private:
    // do not allow copy constructor or assignment
    MsgStateRecv(const MsgStateRecv&);
    operator=(const MsgStateRecv&);
  protected:
    MessageGrp& grp;
    int nbuf; // the number of bytes used in the buffer
    int ibuf; // the current pointer withing the buffer
    int bufsize; // the allocated size of the buffer
    char* buffer; // the data buffer
    char* send_buffer; // the buffer used to send data (includes nbuf)
    int nheader; // nbuf + nheader = the number of bytes in send_buffer to send
    int* nbuf_buffer; // the pointer to the nbuf stored in the buffer

    get_array_void(void*,int);

    virtual void next_buffer() = 0;
  public:
    MsgStateRecv(MessageGrp&);
    virtual ~MsgStateRecv();

    // The buffer size of statein and stateout objects that
    // communicate with each other must match.
    void set_buffer_size(int);

    // I only need to override get(ClassDesc**) but
    // C++ will hide all of the other put's so I must
    // override everything.
    int get(const ClassDesc**);
    int get(char&r);
    int get(int&r);
    int get(float&r);
    int get(double&r);
    int get(char*&);
    int get(int*&);
    int get(float*&);
    int get(double*&);
};

// These allow state to be sent and received and broadcast with
// a message group object.
class StateSend: public MsgStateSend {
  private:
    // do not allow copy constructor or assignment
    StateSend(const StateSend&);
    operator=(const StateSend&);
  private:
    int target_;
  public:
    StateSend(MessageGrp&);
    ~StateSend();
    void target(int);
    void flush();
};

class StateRecv: public MsgStateRecv {
  private:
    // do not allow copy constructor or assignment
    StateRecv(const StateRecv&);
    operator=(const StateRecv&);
  private:
    int source_;
  protected:
    void next_buffer();
  public:
    StateRecv(MessageGrp&);
    void source(int);
};

// Only one node uses a BcastStateSend and the rest must
// use a BcastStateIn.
class BcastStateSend: public MsgStateSend {
  private:
    // do not allow copy constructor or assignment
    BcastStateSend(const BcastStateSend&);
    operator=(const BcastStateSend&);
  public:
    BcastStateSend(MessageGrp&);
    ~BcastStateSend();
    void flush();
};

class BcastStateRecv: public MsgStateRecv {
  private:
    // do not allow copy constructor or assignment
    BcastStateRecv(const BcastStateRecv&);
    operator=(const BcastStateRecv&);
  protected:
    int source_;
    void next_buffer();
  public:
    BcastStateRecv(MessageGrp&, int source = 0);
    void source(int s);
};

#endif
