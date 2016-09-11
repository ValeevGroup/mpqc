
#ifndef SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_
#define SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_

#include "keyval.hpp"

namespace mpqc {
namespace detail {

/** This, together with ForceLink, is used to force code for particular
    classes to be linked into executables. */
template <class A>
class ForceLinkBase {
 public:
  ForceLinkBase() {}
  virtual ~ForceLinkBase() {}
  virtual DescribedClass *create(A) = 0;
};

/** This, together with ForceLinkBase, is used to force code for particular
classes to be linked into executables.  Objects are created from input and
checkpoint files by using class name lookup to find that class's ClassDesc
object.  The ClassDesc object has members that can create the class.
Unfortunately, linking in a library doesn't cause code for the
ClassDesc, and thus the class itself, to be linked.  ForceLink objects are
created in linkage.h files for each library.  The code containing the main
routine for an executable can include these linkage files to force code for
that library's classes to be linked. */
template <class T, class A = const KeyVal& >
class ForceLink : public ForceLinkBase<A> {
 public:
  ForceLink() {}
  virtual ~ForceLink() {}
  DescribedClass *create(A a) { return new T(a); }
};

}  // namespace detail
}  // namespace mpqc

#endif  // SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_
