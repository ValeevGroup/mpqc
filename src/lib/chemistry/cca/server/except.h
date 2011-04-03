#ifndef _except_h
#define _except_h

#include <stdexcept>

class errno_exception: public std::exception {
    std::string what_;
  public:
    ~errno_exception() throw();
    errno_exception(const std::string &wh);
    const char *what() const throw ();
};

#endif
