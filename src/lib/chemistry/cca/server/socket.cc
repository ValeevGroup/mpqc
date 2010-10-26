#ifdef __GNUG__
#pragma implementation
#endif

#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <errno.h>
#include <netinet/in.h>
#include <netdb.h>

#include <iostream>

#include "socket.h"
#include "except.h"

/////////////////////////////////////////////////////////////////////
// TCPSocket

TCPSocket::TCPSocket()
{
  socket_ = 0;
  initialized_ = false;
  bound_ = false;
}

TCPSocket::~TCPSocket()
{
  if (initialized_) close();
}

void
TCPSocket::create()
{
  do { errno = 0; socket_ = socket(PF_INET, SOCK_STREAM, 0); } while (errno == EINTR);
  if (socket_ < 0) throw errno_exception("TCPSocket::create");
  initialized_ = true;
}

void
TCPSocket::bind(u_int16_t port_begin, u_int16_t port_fence)
{
  if (!initialized_)
      throw std::runtime_error("TCPSocket::bind: not initialized");

  struct sockaddr_in ssockaddr;
  ssockaddr.sin_addr.s_addr = INADDR_ANY;
  ssockaddr.sin_family = AF_INET;
  u_int16_t port = port_begin;
  for (port=port_begin; port < port_fence; port++) {
      ssockaddr.sin_port = htons(port);
      if (::bind(socket_, (sockaddr*)&ssockaddr, sizeof(ssockaddr))==0) {
          port_ = port;
          bound_ = true;
          return;
        }
    }
  throw errno_exception("TCPSocket::bind: failed");
}

void
TCPSocket::close()
{
  if (!initialized_)
      throw std::runtime_error("TCPSocket::close: not initialized");

  if (::close(socket_) != 0)
      throw errno_exception("TCPSocket::close");

  initialized_ = false;
  bound_ = false;
}

u_int32_t
TCPSocket::addr()
{
  char hostname[256];
  gethostname(hostname,256);
  struct hostent *hent = gethostbyname(hostname);
  std::cout << "hostname = " << hostname << std::endl;
  u_int32_t add = htonl(*(unsigned long*)hent->h_addr_list[0]);
  std::cout << "TCPSocket::addr() = " << (void*)add << std::endl;
  return add;
}

/////////////////////////////////////////////////////////////////////
// TCPIOSocket

int
TCPIOSocket::read(void *d, int n)
{
  if (!initialized_) throw std::runtime_error("not inited for read");

  int nleft = n;
  while (nleft) {
      int ntransfer = ::read(socket_, d, n);
      if (ntransfer == -1 ) throw errno_exception("socket read");
      nleft -= ntransfer;
    }

  return n;
}

int
TCPIOSocket::write(const void *d, int n)
{
  if (!initialized_) throw std::runtime_error("not inited for write");

  int nleft = n;
  while (nleft) {
      int ntransfer = ::write(socket_, d, n);
      if (ntransfer == -1 ) throw errno_exception("socket write");
      nleft -= ntransfer;
    }

  return n;
}

int
TCPIOSocket::read_string(std::string &s)
{
  int nbytes = 0;
  int size;
  nbytes += read_int(size);
  char *dat = new char[size+1];
  nbytes += read(dat, size);
  dat[size] = '\0';
  s = dat;
  delete[] dat;
  return nbytes;
}

int
TCPIOSocket::write_string(const std::string &s)
{
  int nbytes = 0;
  int size = s.size();
  nbytes += write_int(size);
  nbytes += write(s.c_str(), s.size());
  return nbytes;
}

int
TCPIOSocket::read_int(int &i)
{
  return read(&i, sizeof(int));
}

int
TCPIOSocket::write_int(int i)
{
  return write(&i, sizeof(int));
}

int
TCPIOSocket::read_uint32(u_int32_t &i)
{
  return read(&i, sizeof(u_int32_t));
}

int
TCPIOSocket::write_uint32(u_int32_t i)
{
  return write(&i, sizeof(u_int32_t));
}

/////////////////////////////////////////////////////////////////////
// TCPServerSocket

void
TCPServerSocket::listen(int queue_length)
{
  if (!initialized_)
      throw std::runtime_error("TCPServerSocket::listen: not initialized");

  if (::listen(socket_, queue_length))
      throw errno_exception("TCPServerSocket::listen");

}

/////////////////////////////////////////////////////////////////////
// TCPServerConnection

void
TCPServerConnection::accept(const TCPServerSocket &s)
{
  // this should not be initialized if used with accept
  if (initialized_)
      throw std::runtime_error("TCPServerConnection::accept: already initialized");

  struct sockaddr_in ssockaddr;
  socklen_t len = sizeof(ssockaddr);
  socket_ = ::accept(s.socket_,
                     (sockaddr*)&ssockaddr, &len);
  if (socket_ < 0)
      throw errno_exception("TCPServerConnection::accept");

  initialized_ = true;
}

/////////////////////////////////////////////////////////////////////
// TCPClientConnection

TCPClientConnection::TCPClientConnection()
{
  connected_ = false;
}

void
TCPClientConnection::close()
{
  connected_ = false;
  TCPIOSocket::close();
}

void
TCPClientConnection::connect(const char *remote_hostname,
                             u_int16_t remote_port)
{
  if (!initialized_)
     throw std::runtime_error("TCPClientConnection::connect: not initialized");

  struct hostent * remote_he = gethostbyname(remote_hostname);
  u_int32_t remote_host = htonl(*(unsigned long*)remote_he->h_addr_list[0]);
  connect(remote_host,remote_port);
  connected_ = true;
}

void
TCPClientConnection::connect(u_int32_t remote_host, u_int16_t remote_port)
{
  if (!initialized_)
     throw std::runtime_error("TCPClientConnection::connect: not initialized");

  struct sockaddr_in ssockaddr;
  ssockaddr.sin_addr.s_addr = htonl(remote_host);
  ssockaddr.sin_family = AF_INET;
  ssockaddr.sin_port = htons(remote_port);

  if (::connect(socket_, (struct sockaddr*)&ssockaddr, sizeof(ssockaddr))) {
      throw errno_exception("TCPClientConnection::connect");
    }
  connected_ = true;
}

