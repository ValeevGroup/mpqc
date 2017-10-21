#ifndef MPQC4_SRC_MPQC_UTIL_CORE_BACKTRACE_H_
#define MPQC4_SRC_MPQC_UTIL_CORE_BACKTRACE_H_

#include <string>
#include <vector>

namespace mpqc {
namespace detail {
/**
 * Creates a backtrace of a running program/thread. Example of use:
 * \code
 * void make_omelet(int num_eggs) {
 *   if (num_eggs < 1) {
 *     mpqc::detail::Backtrace bt("breakfast fail:");
 *     throw std::runtime_error(bt.str());
 *   }
 *   stove.on();
 *   // etc.
 * }
 * \endcode
 *
 */
class Backtrace {
 public:
  /**
   * @param prefix will be prepended to each line
   */
  Backtrace(const std::string& prefix = std::string(""));
  Backtrace(const Backtrace&);

  /**
   * @return true is did not get a backtrace
   */
  bool empty() const { return frames_.empty(); }

  /**
   * converts to a string
   * @param nframes_to_skip how many frames to skip
   * @return string representation of Backtrace, with each frame on a separate
   * line, from bottom to top
   */
  std::string str(const size_t nframes_to_skip = 0) const;

 private:
  /// frames_.begin() is the bottom of the stack
  std::vector<std::string> frames_;
  /// prepended to each line
  std::string prefix_;

  /// demangles a symbol
  static std::string __demangle(const std::string& symbol);
};
}  // namespace detail
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_UTIL_CORE_BACKTRACE_H_
