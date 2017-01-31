
#ifndef MPQC4_SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_
#define MPQC4_SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_

#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace detail {

/** This, together with ForceLink, is used to force code for particular
    classes to be linked into executables. */
template <class A>
class ForceLinkBase {
 public:
  ForceLinkBase() {}
  virtual ~ForceLinkBase() {}
  virtual DescribedClass* create(A) = 0;
};

/**
 * \brief This, together with ForceLinkBase, is used to force code for
 * particular
 *        classes to be linked into executables.
 *
 * \note Objects of classes derived from DescribedClass can be created from
 * KeyVal input
 *       and/or checkpoint files by using class name lookup to find that class's
 * entry
 *       in the constructor registry. Unfortunately, linking in a library that
 * defines
 *       such a class doesn't cause the code for the class to be linked. To
 * force linking
 *       of class T derived from DescribedClass create an object of class
 * ForceLink<T>
 *       in the main executable. This is most conveniently done by having every
 * library
 *       define such objects in its \c linkage.h file; the main executable
 *       then includes all \c linkage.h files from the libraries that it needs.
 * \note Only a declaration of class T is needed to define an object of type
 * ForceLink<T>;
 *       this is key to avoiding instantiation of class T in the executable.
 *       Nevertheless, it is necessary to define ForceLink<T>::create in the
 * library somewhere,
 *       usually in a file called \c linkage.cpp ; this file usually also
 * contains
 *       constructor registration macros.
 * \sa MPQC_FORCELINK_KEYVAL_CTOR
 */
template <class T, class A = const KeyVal&>
class ForceLink : public ForceLinkBase<A> {
 public:
  ForceLink() {}
  virtual ~ForceLink() {}
  DescribedClass* create(A a);
};

}  // namespace detail
}  // namespace mpqc

/// Streamlines the definition of ForceLink<T>::create methods

/// For every object of ForceLink<T> type defined in \c linkage.h
/// add MPQC_FORCELINK_KEYVAL_CTOR(T) to a .cpp file.
/// \note Since for DescribedClass-based classes one needs to
///       register the KeyVal ctor AND force their linking
///       into the executable, usually this macro should not be used
///       directly, instead use MPQC_CLASS_EXPORT2(Key,T) or
///       MPQC_CLASS_EXPORT(T)
/// \sa MPQC_CLASS_EXPORT
/// \sa MPQC_CLASS_EXPORT2
#define MPQC_FORCELINK_KEYVAL_CTOR(...)                          \
  namespace mpqc {                                               \
  namespace detail {                                             \
  template <>                                                    \
  DescribedClass* ForceLink<__VA_ARGS__, const KeyVal&>::create( \
      const KeyVal& kv) {                                        \
    return new __VA_ARGS__(kv);                                  \
  }                                                              \
  }                                                              \
  }
/**/

/// For a class \c T derived from DescribedClass
/// MPQC_CLASS_EXPORT(T), creates a Boost.Serialization GUID for the class,
/// registers the KeyVal ctor of \c T (see MPQC_CLASS_EXPORT_KEY(T)),
/// and
/// allows to force linking of \c T into the executable by defining
/// a ForceLink<T> object in it (see MPQC_FORCELINK_KEYVAL_CTOR(T)).
/// Place MPQC_CLASS_EXPORT(T) to a source file in global scope.
/// \note To avoid slowing down compilation by coupling the instantiation
///       of multiple classes, use only one MPQC_CLASS_EXPORT(T) per .cpp file
#define MPQC_CLASS_EXPORT(...)       \
  MPQC_CLASS_EXPORT_KEY(__VA_ARGS__) \
  MPQC_FORCELINK_KEYVAL_CTOR(__VA_ARGS__)

/// For a class \c T derived from DescribedClass
/// MPQC_CLASS_EXPORT2(Key,T) creates a Boost.Serialization
/// GUID for the class using key \c Key
/// registers the KeyVal ctor of \c T (see MPQC_CLASS_EXPORT_KEY2(Key,T)),
/// and
/// allows to force linking of \c T into the executable by defining
/// a ForceLink<T> object in it (see MPQC_FORCELINK_KEYVAL_CTOR(T)).
/// Place MPQC_CLASS_EXPORT2(T) to a source file in global scope.
/// \note To avoid slowing down compilation by coupling the instantiation
///       of multiple classes, use only one MPQC_CLASS_EXPORT2(T) per .cpp file
#define MPQC_CLASS_EXPORT2(K, ...)       \
  MPQC_CLASS_EXPORT_KEY2(K, __VA_ARGS__) \
  MPQC_FORCELINK_KEYVAL_CTOR(__VA_ARGS__)

#endif  // MPQC4_SRC_MPQC_UTIL_KEYVAL_FORCELINK_H_
