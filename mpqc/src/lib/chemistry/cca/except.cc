#ifdef __GNUG__
#pragma implementation
#endif

#include <errno.h>

#include "except.h"

errno_exception::~errno_exception() throw()
{
}

errno_exception::errno_exception(const std::string &wh)
{
  what_ = wh + ": " + strerror(errno);
}

const char *errno_exception::what() const throw ()
{
  return what_.c_str();
}
