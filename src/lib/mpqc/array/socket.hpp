#ifndef MPQC_ARRAY_SOCKET_HPP
#define MPQC_ARRAY_SOCKET_HPP

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#include <limits.h>
#include <ifaddrs.h>

namespace mpqc {
namespace detail {

    struct ArraySocket {

        typedef sockaddr Address;

        static void assert_(bool cond, const std::string &msg) {
            if (!cond) {
                perror(msg.c_str());
                throw std::runtime_error(msg + " failed");
            }
        }

        ArraySocket() {
            this->fd_ = -1;
        }

        void start() {

            // sockaddr addr;
            // {
            //     struct ifaddrs *ifs;
            //     if (getifaddrs(&ifs) == -1) {
            //         perror("getifaddrs");
            //         exit(EXIT_FAILURE);
            //     }
            //     /* Walk through linked list, maintaining head pointer so we
            //        can free list later */
            //     for (auto ifa = ifs; ifa != NULL; ifa = ifa->ifa_next) {
            //         if (ifa->ifa_addr == NULL) continue;
            //         if (ifa->ifa_addr->sa_family == AF_INET) {
            //             struct sockaddr_in *sin = (struct sockaddr_in*)ifa->ifa_addr;
            //             printf("Interface %s, addr=%s\n",
            //                    ifa->ifa_name, inet_ntoa(sin->sin_addr));
            //         }
            //     }
            //     freeifaddrs(ifs);
            // }

            struct in_addr in;
            {
                int err;
                char hostname[512];
                struct hostent *h;
                assert_(::gethostname(hostname, sizeof(hostname)) > -1, "gethostname");
                assert_((h = ::gethostbyname(hostname)), "gethostbyname");
                if (h->h_addrtype == AF_INET) {
                    in = *(in_addr*)h->h_addr;
                    // //err = inet_pton(AF_INET, h->h_addr, &in);
                    // // if (err <= 0) perror("inet_pton");
                    // printf("using address %s(%i)\n", inet_ntoa(in), in.s_addr);
                }
            }

            struct sockaddr_in sin;
            sin.sin_family = AF_INET;
            sin.sin_port = 0;
            sin.sin_addr.s_addr = in.s_addr;
            memset(sin.sin_zero, '\0', sizeof(sin.sin_zero));

            int sock;
            assert_((sock = ::socket(PF_INET, SOCK_STREAM, 0)) > -1, "socket");
            assert_(::bind(sock, (struct sockaddr*)&sin, sizeof(sin))  > -1, "bind");
            assert_(::listen(sock, 10) > -1, "listen");

            this->fd_ = sock;

            {
                Address addr = this->address();
                printf("ArrayServer running on %s:%i\n",
                       inet_ntoa(((sockaddr_in*)&addr)->sin_addr),
                       ((sockaddr_in*)&addr)->sin_port);
            }

        }

        Address address() const {
            struct sockaddr addr;
            socklen_t len = sizeof(addr);
            assert_(getsockname(this->fd_, &addr, &len) > -1, "getsockname");
            return addr;
        }

        template<typename T>
        void wait(T *data) const {
            struct sockaddr_storage addr;
            socklen_t addrlen = sizeof(addr);
            int fd, bytes;
            assert_((fd = ::accept(this->fd_, (sockaddr*)&addr, &addrlen)) > -1,
                    "accept");
            assert_((bytes = ::recv(fd, data, sizeof(T), 0)) > -1, "recv");
            assert_(::close(fd) > -1, "close");
            if (bytes != sizeof(T)) {
                throw std::runtime_error("wrong number of bytes received");
            }
        }

        template<typename T>
        static void send(const T *data, struct sockaddr addr) {
            int fd;
            assert_((fd = ::socket(AF_INET, SOCK_STREAM, 0)) > -1, "socket");
            assert_(::connect(fd, &addr, sizeof(addr)) > -1, "connect");
            assert_(::send(fd, data, sizeof(T), 0) > -1, "send");
            assert_(::close(fd) > -1, "close");
        }

    private:
        int fd_;

    };

} // namespace detail
} // namespace mpqc

#endif // MPQC_ARRAY_SOCKET_HPP
