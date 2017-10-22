#include "backtrace.h"

#include <iterator>
#include <sstream>

#include "mpqc/mpqc_config.h"
#if defined(HAVE_LIBUNWIND)
#define UNW_LOCAL_ONLY
#include <libunwind.h>
#elif defined(HAVE_BACKTRACE)
#include <execinfo.h>
#endif

#ifdef HAVE_CXA_DEMANGLE
#include <cxxabi.h>
#endif

#include "formio.h"

namespace mpqc {
namespace detail {
Backtrace::Backtrace(const std::string &prefix) : prefix_(prefix) {
#ifdef HAVE_LIBUNWIND
  {
    unw_cursor_t cursor;
    unw_context_t uc;
    unw_word_t ip, sp, offp;
    int frame = 0;

    unw_getcontext(&uc);
    unw_init_local(&cursor, &uc);
    while (unw_step(&cursor) > 0) {
      unw_get_reg(&cursor, UNW_REG_IP, &ip);
      unw_get_reg(&cursor, UNW_REG_SP, &sp);
      char name[1024];
      unw_get_proc_name(&cursor, name, 1024, &offp);
      std::ostringstream oss;
      oss << prefix_ << "frame " << frame << ": "
          << mpqc::printf("ip = 0x%lx sp = 0x%lx ", (long)ip, (long)sp)
          << " symbol = " << __demangle(name);
      frames_.push_back(oss.str());
      ++frame;
    }
  }
#elif defined(HAVE_BACKTRACE)  // !HAVE_LIBUNWIND
  void *stack_addrs[1024];
  const int naddrs = backtrace(stack_addrs, 1024);
  char **frame_symbols = backtrace_symbols(stack_addrs, naddrs);
  // starting @ 1 to skip this function
  for (int i = 1; i < naddrs; ++i) {
    // extract (mangled) function name
    std::string mangled_function_name;
    {
      std::istringstream iss(std::string(frame_symbols[i]),
                             std::istringstream::in);
      std::string frame, file, address;
      iss >> frame >> file >> address >> mangled_function_name;
    }

    std::ostringstream oss;
    oss << prefix_ << "frame " << i << ": return address = " << stack_addrs[i]
        << std::endl
        << "  symbol = " << __demangle(mangled_function_name);
    frames_.push_back(oss.str());
  }
  free(frame_symbols);
#else                          // !HAVE_LIBUNWIND && !HAVE_BACKTRACE
#if defined(SIMPLE_STACK)
  int bottom = 0x1234;
  void **topstack = (void **)0xffffffffL;
  void **botstack = (void **)0x70000000L;
  // signal handlers can put weird things in the return address slot,
  // so it is usually best to keep toptext large.
  void **toptext = (void **)0xffffffffL;
  void **bottext = (void **)0x00010000L;
#endif  // SIMPLE_STACK

#if (defined(linux) && defined(i386))
  topstack = (void **)0xc0000000;
  botstack = (void **)0xb0000000;
#endif
#if (defined(__OSF1__) && defined(i860))
  topstack = (void **)0x80000000;
  botstack = (void **)0x70000000;
#endif

#if defined(SIMPLE_STACK)
  // This will go through the stack assuming a simple linked list
  // of pointers to the previous frame followed by the return address.
  // It trys to be careful and avoid creating new execptions, but there
  // are no guarantees.
  void **stack = (void **)&bottom;

  void **frame_pointer = (void **)stack[3];
  while (frame_pointer >= botstack && frame_pointer < topstack &&
         frame_pointer[1] >= bottext && frame_pointer[1] < toptext) {
    std::ostringstream oss;
    oss << prefix_ << "frame: " << (void *)frame_pointer;
    oss << "  retaddr: " << frame_pointer[1];
    frames_.push_back(oss.str());

    frame_pointer = (void **)*frame_pointer;
  }
#endif  // SIMPLE_STACK
#endif  // HAVE_BACKTRACE
}

Backtrace::Backtrace(const Backtrace &other)
    : frames_(other.frames_), prefix_(other.prefix_) {}

std::string Backtrace::str(size_t nframes_to_skip) const {
  std::ostringstream oss;
  std::copy(frames_.begin() + nframes_to_skip, frames_.end(),
            std::ostream_iterator<std::string>(oss, "\n"));
  return oss.str();
}

std::string Backtrace::__demangle(const std::string &symbol) {
  std::string dsymbol;
#ifdef HAVE_CXA_DEMANGLE
  {
    int status;
    char *dsymbol_char = abi::__cxa_demangle(symbol.c_str(), 0, 0, &status);
    if (status == 0) {  // success
      dsymbol = dsymbol_char;
      free(dsymbol_char);
    } else  // fail
      dsymbol = symbol;
  }
#else
  dsymbol = symbol;
#endif
  return dsymbol;
}

}  // namespace detail
}  // namespace mpqc
