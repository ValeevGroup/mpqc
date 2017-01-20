#ifndef SRC_MPQC_UTIL_MISC_TASK_H_
#define SRC_MPQC_UTIL_MISC_TASK_H_

namespace mpqc {

/// this is the base for all tasks that the MPQC can execute
/// @internal analogous to Runnable in MPQC3
class Task : public DescribedClass {
  /// Executes an action as specified in the derived class.
  virtual void execute() = 0;
};

}  // namespace mpqc

#endif  // SRC_MPQC_UTIL_MISC_TASK_H_
