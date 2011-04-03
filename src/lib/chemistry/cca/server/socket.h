#ifndef _socket_h
#define _socket_h

#include <unistd.h>
#include <sys/types.h>

#include <string>

class TCPSocket {
  protected:
    TCPSocket();
    virtual ~TCPSocket();

    int socket_;
    int port_;
    bool initialized_;
    bool bound_;

  public:
    void create();
    void bind(u_int16_t port_begin, u_int16_t port_fence);
    int port() { return port_; }
    u_int32_t addr();

    bool initialized() { return initialized_; }
    bool bound() { return bound_; }

    virtual void close();
};

class TCPIOSocket: public TCPSocket {
  public:
    int read(void *d, int n);
    int write(const void *d, int n);
    int read_int(int *d, int n) {return read((void*)d,n*sizeof(int));}
    int write_int(const int *d, int n) {return write((void*)d,n*sizeof(int));}
    int read_string(std::string &);
    int write_string(const std::string &);
    int read_int(int &);
    int write_int(int);
    int read_uint32(u_int32_t &);
    int write_uint32(u_int32_t);
};

class TCPServerSocket: public TCPSocket {
    friend class TCPServerConnection;
  public:
    void listen(int queue_length = 8);
};

class TCPServerConnection: public TCPIOSocket {
  public:
    void accept(const TCPServerSocket &);
};

class TCPClientConnection: public TCPIOSocket {
    bool connected_;
  public:
    TCPClientConnection();
    void close();
    bool connected() { return connected_; }
    void connect(const char *remote_hostname, u_int16_t remote_port);
    void connect(u_int32_t remote_host, u_int16_t remote_port);
};

#endif
