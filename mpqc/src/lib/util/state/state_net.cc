
/* state_net.cc -- implementation of socket StateOut classes
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
#pragma implementation
#endif

#include <iostream.h>
#include <string.h>
#include <util/unix/cct_cprot.h>
#include <util/class/class.h>
#include <util/state/state_net.h>

/*
 * include <netdb.h>  for struct hostent
 * include <sys/socket.h>  for AF_INET, etc
 */

#include <netdb.h>

extern "C" {
#include <sys/types.h>
#include <sys/socket.h>
}

///////////////////////////////////////////////////////////////////

Socket::Socket() :
  fd_(-1)
{
}

Socket::~Socket()
{
  close_socket();
}

int Socket::close_socket()
{
  int type,len=sizeof(int);
  int fd=fd_; fd_=-1;
  if(!getsockopt(fd,SOL_SOCKET,SO_TYPE,(char*)&type,&len)) return close(fd);
  return 0;
  }

///////////////////////////////////////////////////////////////////
/*
 * the following two functions are utility routines for network i/o
 *
 * they are copied pretty much verbatim from:
 *
 *  UNIX Network Programming
 *  W. Richard Stevens
 *  Prentice Hall (Englewood Cliffs, NJ) 1990
 *  
 *  Chapter 6: Berkeley Sockets
 *
 */

/*
 * read "n" bytes from a descriptor.
 */

StreamSocket::StreamSocket() :
  Socket()
{
}

StreamSocket::~StreamSocket()
{
}

int StreamSocket::readn(void *ptr, int nbytes)
{
  int nleft=nbytes;
  while(nleft > 0) {
    int nread = read(fd_,(char*)ptr,nleft);
    if (nread < 0) return(nread);   // error, return < 0
    else if (nread==0) break;       // EOF

    nleft -= nread;
#ifdef __GNUC__
    ptr += nread;
#else
    char *ppp=(char*)ptr;
    ppp += nread;
    ptr=ppp;
#endif
    }

  return(nbytes-nleft);     // return >= 0
  }

/*
 * Write "n" bytes from a descriptor.
 */
int StreamSocket::writen(void *ptr, int nbytes)
{
  int nleft=nbytes;
  while(nleft > 0) {
    int nwritten = write(fd_,(char*)ptr,nleft);
    if (nwritten <= 0) return(nwritten);   // error, return <= 0

    nleft -= nwritten;
#ifdef __GNUC__
    ptr += nwritten;
#else
    char *ppp=(char*)ptr;
    ppp += nwritten;
    ptr=ppp;
#endif
    }

  return(nbytes-nleft);     // return >= 0
  }

///////////////////////////////////////////////////////////////////

DGramSocket::DGramSocket() :
  Socket()
{
}

DGramSocket::~DGramSocket()
{
}

// send a datagram of n bytes
int DGramSocket::sendn(void *p, int n, struct sockaddr *raddr, int addrlen)
{
  int nsent=sendto(fd_,(char*)p,n,0,raddr,addrlen);
  return nsent;
  }

// receive a datagram of n bytes
int DGramSocket::recvn(void *p, int n, struct sockaddr *raddr, int& addrlen)
{
  int nrcvd=recvfrom(fd_,(char*)p,n,0,raddr,&addrlen);
  return nrcvd;
  }

///////////////////////////////////////////////////////////////////

BSDSocket::BSDSocket() :
  Socket()
{
}

BSDSocket::~BSDSocket()
{
}

// bind the local address to port "port"
int BSDSocket::bind_socket(int port)
{
  memset((char*)&laddr_,'\0',sizeof(laddr_));

  laddr_.sin_family = AF_INET;
  laddr_.sin_addr.s_addr = htonl(INADDR_ANY);
  laddr_.sin_port = htons(port);

  if(bind(fd_,(struct sockaddr *) &laddr_, sizeof(laddr_)) < 0) {
    err_ret("BSDSocket::bind_socket():  could not bind local address");
    return -1;
    }
  return 0;
  }

// make connection to host "host" on port "port"
int BSDSocket::connect_socket(const char *host,int port)
{
  memset((char*)&raddr_,'\0',sizeof(raddr_));
  raddrlen=sizeof(raddr_);

  struct hostent *he = gethostbyname(host);
  if(he==0) {
    err_ret("BSDSocket::connect_socket(): unknown host %s",host);
    return -1;
    }

  raddr_.sin_family=he->h_addrtype;
  raddr_.sin_addr.s_addr = ((struct in_addr *) he->h_addr_list[0])->s_addr;
  raddr_.sin_port = htons(port);

  int type,len=sizeof(int);
  if(getsockopt(fd_,SOL_SOCKET,SO_TYPE,(char*)&type,&len)) {
    err_ret("BSDSocket::connect_socket(): cannot determine type of socket");
    return -1;
    }
    
 // only call connect() for TCP sockets, let recvfrom do the work for UDP
  if(type==SOCK_STREAM) {
    if(connect(fd_,(struct sockaddr *)&raddr_, raddrlen) < 0) {
      err_ret("BSDSocket::connect_socket(): could not connect to host %s",host);
      return -1;
      }
    }
  return 0;
  }

///////////////////////////////////////////////////////////////////

BSD_TCPSocket::BSD_TCPSocket()
{
  if((fd_ = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    err_sys("BSD_TCPSocket(): could not open stream socket");
  }

BSD_TCPSocket::~BSD_TCPSocket()
{
}

// bind the local address, and then do a listen call to tell the system
// that this is a connection-oriented server
int BSD_TCPSocket::bind_socket(int port)
{
  if(BSDSocket::bind_socket(port)<0) return -1;
  return listen(fd_,5);
  }

// wait for connection from clients, returns file descriptor for socket to
//   client
// this is called Listen instead of Accept since the TLI equivalent is
//   t_listen()
int BSD_TCPSocket::listen_socket()
{
  raddrlen=sizeof(raddr_);
  return accept(fd_,(struct sockaddr *) &raddr_, &raddrlen);
  }

int BSD_TCPSocket::listen_socket(BSD_TCPSocket& clnt)
{
  clnt.raddrlen=sizeof(clnt.raddr_);
  clnt.fd_ = accept(fd_,(struct sockaddr *) &clnt.raddr_, &clnt.raddrlen);
  return clnt.fd_;
  }

///////////////////////////////////////////////////////////////////

BSD_UDPSocket::BSD_UDPSocket()
{
  if((fd_ = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    err_sys("BSD_UDPSocket(): could not open datagram socket");
  }

BSD_UDPSocket::~BSD_UDPSocket()
{
}

///////////////////////////////////////////////////////////////////

StateOutBSD_TCP::StateOutBSD_TCP()
{
}

StateOutBSD_TCP::~StateOutBSD_TCP()
{
}

int StateOutBSD_TCP::put_array_void(void* p, int n)
{
  if(writen(p,n)!=n) {
    err_sys("StateOutBSD_TCP::put_array_void(void* p, int n): "
            "trouble writing data to socket");
    }
  return n;
  }

StateInBSD_TCP::StateInBSD_TCP()
{
}

StateInBSD_TCP::~StateInBSD_TCP()
{
}

int StateInBSD_TCP::get_array_void(void* p, int n)
{
  if(readn(p,n)!=n) {
    err_sys("StateInBSD_TCP::get_array_void(void* p, int n): "
            "trouble reading data from socket");
    }
  return n;
  }

///////////////////////////////////////////////////////////////////

StateIOBSD_TCP::StateIOBSD_TCP()
{
}

StateIOBSD_TCP::~StateIOBSD_TCP()
{
}

///////////////////////////////////////////////////////////////////

StateOutBSD_UDP::StateOutBSD_UDP()
{
}

StateOutBSD_UDP::~StateOutBSD_UDP()
{
}

int StateOutBSD_UDP::put_array_void(void* p, int n)
{
  if(sendn(p,n,(struct sockaddr*)&raddr_,raddrlen)!=n) {
    err_sys("StateOutBSD_UDP::put_array_void(void* p, int n): "
            "trouble writing data to socket");
    }
  return n;
  }

StateInBSD_UDP::StateInBSD_UDP()
{
}

StateInBSD_UDP::~StateInBSD_UDP()
{
}

int StateInBSD_UDP::get_array_void(void* p, int n)
{
  if(recvn(p,n,(struct sockaddr *)&raddr_,raddrlen)!=n) {
    err_sys("StateInBSD_UDP::get_array_void(void* p, int n): "
            "trouble reading data from socket");
    }
  return n;
  }

///////////////////////////////////////////////////////////////////

StateIOBSD_UDP::StateIOBSD_UDP()
{
}

StateIOBSD_UDP::~StateIOBSD_UDP()
{
}

///////////////////////////////////////////////////////////////////

UnixSocket::UnixSocket()
{
}

UnixSocket::~UnixSocket()
{
  unlink_socket();
}

int UnixSocket::unlink_socket()
{
  Socket::close_socket();
  return unlink(laddr_.sun_path);
  }

// bind the local address to the socket "path"
int UnixSocket::bind_socket(const char *path)
{
  memset((char*)&laddr_,'\0',sizeof(laddr_));

  laddr_.sun_family = AF_UNIX;
  if(path)
    strncpy(laddr_.sun_path,path,sizeof(laddr_.sun_path));
  else {
    strcpy(laddr_.sun_path,"/tmp/dg.XXXXXX");
    mktemp(laddr_.sun_path);
    }
  int len=strlen(laddr_.sun_path)+sizeof(laddr_.sun_family);

  if(bind(fd_,(struct sockaddr *) &laddr_, len) < 0) {
    err_ret("UnixSocket::bind_socket():  could not bind local address");
    return -1;
    }
  return 0;
  }

// make connection to socket "path"
int UnixSocket::connect_socket(const char *path)
{
  memset((char*)&raddr_,'\0',sizeof(raddr_));

  raddr_.sun_family=AF_UNIX;
  strncpy(raddr_.sun_path,path,sizeof(raddr_.sun_path));
  raddrlen=strlen(raddr_.sun_path)+sizeof(raddr_.sun_family);

  int type=0,len=sizeof(int);
  if(getsockopt(fd_,SOL_SOCKET,SO_TYPE,(char*)&type,&len)) {
    err_ret("UnixSocket::connect_socket(): cannot determine type of socket");
    return -1;
    }
    
 // only call connect() for stream sockets, let recvfrom do the work for dgram
  if(type==SOCK_STREAM) {
    if(connect(fd_,(struct sockaddr *)&raddr_, raddrlen) < 0) {
      err_ret("UnixSocket::connect_socket(): could not connect to socket %s",path);
      return -1;
      }
    }
  return 0;
  }

///////////////////////////////////////////////////////////////////

UnixStreamSocket::UnixStreamSocket()
{
  if((fd_ = socket(AF_UNIX, SOCK_STREAM, 0)) < 0)
    err_sys("UnixStreamSocket(): could not open stream socket");
  }

UnixStreamSocket::~UnixStreamSocket()
{
}

// bind the local address, and then do a listen call to tell the system
// that this is a connection-oriented server
int UnixStreamSocket::bind_socket(const char *host)
{
  if(UnixSocket::bind_socket(host)<0) return -1;
  return listen(fd_,5);
  }

// wait for connection from clients, returns file descriptor for socket to
//   client
// this is called Listen instead of Accept since the TLI equivalent is
//   t_listen()
int UnixStreamSocket::listen_socket()
{
  raddrlen=sizeof(raddr_);
  return accept(fd_,(struct sockaddr *) &raddr_, &raddrlen);
  }

int UnixStreamSocket::listen_socket(UnixStreamSocket& clnt)
{
  clnt.raddrlen=sizeof(clnt.raddr_);
  clnt.fd_ = accept(fd_,(struct sockaddr *) &clnt.raddr_, &clnt.raddrlen);
  return clnt.fd_;
  }

///////////////////////////////////////////////////////////////////

UnixDGramSocket::UnixDGramSocket()
{
  if((fd_ = socket(AF_UNIX, SOCK_DGRAM, 0)) < 0)
    err_sys("UnixDGramSocket(): could not open datagram socket");
  }

UnixDGramSocket::~UnixDGramSocket()
{
}

///////////////////////////////////////////////////////////////////

StateOutUnixStream::StateOutUnixStream()
{
}

StateOutUnixStream::~StateOutUnixStream()
{
}

int StateOutUnixStream::put_array_void(void* p, int n)
{
  if(writen(p,n)!=n) {
    err_sys("StateOutUnixStream::put_array_void(void* p, int n): "
            "trouble writing data to socket");
    }
  return n;
  }

StateInUnixStream::StateInUnixStream()
{
}

StateInUnixStream::~StateInUnixStream()
{
}

int StateInUnixStream::get_array_void(void* p, int n)
{
  if(readn(p,n)!=n) {
    err_sys("StateInUnixStream::get_array_void(void* p, int n): "
            "trouble reading data from socket");
    }
  return n;
  }

///////////////////////////////////////////////////////////////////

StateIOUnixStream::StateIOUnixStream()
{
}

StateIOUnixStream::~StateIOUnixStream()
{
}

///////////////////////////////////////////////////////////////////

StateOutUnixDGram::StateOutUnixDGram()
{
}

StateOutUnixDGram::~StateOutUnixDGram()
{
}

int StateOutUnixDGram::put_array_void(void* p, int n)
{
  if(sendn(p,n,(struct sockaddr*)&raddr_,raddrlen)!=n) {
    err_sys("StateOutUnixDGram::put_array_void(void* p, int n): "
            "trouble writing data to socket");
    }
  return n;
  }

StateInUnixDGram::StateInUnixDGram()
{
}

StateInUnixDGram::~StateInUnixDGram()
{
}

int StateInUnixDGram::get_array_void(void* p, int n)
{
  if(recvn(p,n,(struct sockaddr *)&raddr_,raddrlen)!=n) {
    err_sys("StateInUnixDGram::get_array_void(void* p, int n): "
            "trouble reading data from socket");
    }
  return n;
  }

StateIOUnixDGram::StateIOUnixDGram()
{
}

StateIOUnixDGram::~StateIOUnixDGram()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
