
#ifndef _libqc_state_net_h
#define _libqc_state_net_h

/* state_net.h -- declarations for the socket StateOut classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      October, 1992
 */

#ifdef __GNUC__
#pragma interface
#endif


#include <util/state/state.h>

/*
 * include <unistd.h>     for close()
 * include <sys/types.h>   just in case rpc does not include this
 * include <netinet/in.h>  for struct sockaddr and other IP things
 * include <sys/un.h>      for struct sockaddr_un
 */

#include <unistd.h>
#include <sys/types.h>
#include <netinet/in.h>
#if !defined(I860)
#include <sys/un.h>
#endif

/*
 * include <malloc.h> because the ATT headers on the sun are brain dead
 */
#if !defined(__GNUC__) && defined(SUN)
#include <malloc.h>
#endif

namespace sc {

// ////////////////////////////////////////////////////////////

class Socket {
  protected:
    int fd_;
    Socket();
  public:
    virtual ~Socket();
    virtual int close_socket();
    inline const int fd() const { return fd_; }
  };

// ////////////////////////////////////////////////////////////

class StreamSocket : virtual public Socket {
  protected:
    int writen(void*,int);
    int readn(void*,int);
    StreamSocket();
  public:
    ~StreamSocket();
  };

class DGramSocket : virtual public Socket {
  protected:
    int sendn(void*,int,struct sockaddr*,int);
    int recvn(void*,int,struct sockaddr*,int&);
    DGramSocket();
  public:
    ~DGramSocket();
  };

// ////////////////////////////////////////////////////////////

class BSDSocket : virtual public Socket {
  protected:
    struct sockaddr_in laddr_;
    struct sockaddr_in raddr_;
    int raddrlen;
    BSDSocket();
  public:
    ~BSDSocket();
    virtual int bind_socket(int);
    virtual int connect_socket(const char*,int);
  };

class BSD_TCPSocket: public BSDSocket, public StreamSocket {
  protected:
    BSD_TCPSocket();
  public:
    ~BSD_TCPSocket();

    int bind_socket(int);
    int listen_socket();
    int listen_socket(BSD_TCPSocket&);
  };

class BSD_UDPSocket: public BSDSocket, public DGramSocket {
  protected:
    BSD_UDPSocket();
  public:
    ~BSD_UDPSocket();
  };

// ////////////////////////////////////////////////////////////

class StateOutBSD_TCP : virtual public BSD_TCPSocket, public StateOutXDR {
  protected:
    int put_array_void(const void*,int);
  public:
    StateOutBSD_TCP();
    ~StateOutBSD_TCP();
  };

class StateInBSD_TCP : virtual public BSD_TCPSocket, public StateInXDR {
  protected:
    int get_array_void(void*,int);
  public:
    StateInBSD_TCP();
    ~StateInBSD_TCP();
  };

class StateIOBSD_TCP : public StateInBSD_TCP, public StateOutBSD_TCP {
  public:
    StateIOBSD_TCP();
    ~StateIOBSD_TCP();
  };

// ////////////////////////////////////////////////////////////

class StateOutBSD_UDP : virtual public BSD_UDPSocket, public StateOutXDR {
  protected:
    int put_array_void(const void*,int);
  public:
    StateOutBSD_UDP();
    ~StateOutBSD_UDP();
  };

class StateInBSD_UDP : virtual public BSD_UDPSocket, public StateInXDR {
  protected:
    int get_array_void(void*,int);
  public:
    StateInBSD_UDP();
    ~StateInBSD_UDP();
  };

class StateIOBSD_UDP : public StateInBSD_UDP, public StateOutBSD_UDP {
  public:
    StateIOBSD_UDP();
    ~StateIOBSD_UDP();
  };

// ////////////////////////////////////////////////////////////

class UnixSocket : virtual public Socket {
  protected:
    struct sockaddr_un laddr_;
    struct sockaddr_un raddr_;
    int raddrlen;
    UnixSocket();
  public:
    virtual ~UnixSocket();
    int unlink_socket();
    virtual int bind_socket(const char*);
    virtual int connect_socket(const char*);
  };

class UnixStreamSocket: public UnixSocket, public StreamSocket {
  protected:
    UnixStreamSocket();
  public:
    ~UnixStreamSocket();

    int bind_socket(const char*);
    int listen_socket();
    int listen_socket(UnixStreamSocket&);
  };

class UnixDGramSocket: public UnixSocket, public DGramSocket {
  protected:
    UnixDGramSocket();
  public:
    ~UnixDGramSocket();
  };

// ////////////////////////////////////////////////////////////

class StateOutUnixStream : virtual public UnixStreamSocket, public StateOut {
  protected:
    int put_array_void(const void*,int);
  public:
    StateOutUnixStream();
    ~StateOutUnixStream();
  };

class StateInUnixStream : virtual public UnixStreamSocket, public StateIn {
  protected:
    int get_array_void(void*,int);
  public:
    StateInUnixStream();
    ~StateInUnixStream();
  };

class StateIOUnixStream: public StateInUnixStream, public StateOutUnixStream {
  public:
    StateIOUnixStream();
    ~StateIOUnixStream();
  };

// ////////////////////////////////////////////////////////////

class StateOutUnixDGram : virtual public UnixDGramSocket, public StateOut {
  protected:
    int put_array_void(const void*,int);
  public:
    StateOutUnixDGram();
    ~StateOutUnixDGram();
  };

class StateInUnixDGram : virtual public UnixDGramSocket, public StateIn {
  protected:
    int get_array_void(void*,int);
  public:
    StateInUnixDGram();
    ~StateInUnixDGram();
  };

class StateIOUnixDGram : public StateInUnixDGram, public StateOutUnixDGram {
  public:
    StateIOUnixDGram();
    ~StateIOUnixDGram();
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
