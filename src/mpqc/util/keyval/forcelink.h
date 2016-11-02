
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

/**
 * \brief This, together with ForceLinkBase, is used to force code for particular
 *        classes to be linked into executables.
 *
 * \note Objects of classes derived from DescribedClass can be created from KeyVal input
 *       and/or checkpoint files by using class name lookup to find that class's entry
 *       in the constructor registry. Unfortunately, linking in a library that defines
 *       such a class doesn't cause the code for the class to be linked. To force linking
 *       of class T derived from DescribedClass create an object of class ForceLink<T>
 *       in the main executable. This is most conveniently done by having every library
 *       define such objects in its \c linkage.h file; the main executable
 *       then includes all \c linkage.h files from the libraries that it needs.
 * \note Only a declaration of class T is needed to define an object of type ForceLink<T>;
 *       this is key to avoiding instantiation of class T in the executable.
 *       Nevertheless, it is necessary to define ForceLink<T>::create in the library somewhere,
 *       usually in a file called \c linkage.cpp ; this file usually also contains
 *       constructor registration macros.
 * \sa MPQC_FORCELINK_KEYVAL_CTOR
 */
template <class T, class A = const KeyVal& >
class ForceLink : public ForceLinkBase<A> {
 public:
  ForceLink() {}
  virtual ~ForceLink() {}
  DescribedClass *create(A a);
};

}  // namespace detail
}  // namespace mpqc


/// Streamlines the definition of ForceLink<T>::create methods

/// For every object of ForceLink<T> type defined in \c linkage.h
/// add MPQC_FORCELINK_KEYVAL_CTOR(T) to \c linkage.cpp file.
#define MPQC_FORCELINK_KEYVAL_CTOR(...)                              \
  template <>                                                        \
  mpqc::DescribedClass*                                              \
  mpqc::detail::ForceLink<__VA_ARGS__, const mpqc::KeyVal&>::create( \
      const mpqc::KeyVal& kv) {                                      \
    return new __VA_ARGS__(kv);                                      \
  }

#endif  // SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_
