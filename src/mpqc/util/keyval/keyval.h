/*
 * keyval.h
 *
 *  Created on: Apr 18, 2016
 *      Author: evaleev
 */

#ifndef MPQC4_SRC_MPQC_UTIL_KEYVAL_KEYVAL_H_
#define MPQC4_SRC_MPQC_UTIL_KEYVAL_KEYVAL_H_

#include <array>
#include <cassert>
#include <list>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/serialization/export.hpp>

#include "mpqc/util/meta/predicates.h"
#include "mpqc/util/meta/stream.h"

// serialize all pointers as void*
// NB XCode 7.3.1 (7D1014) libc++ char stream does not properly deserialize
// void*, use size_t instead
namespace boost {
namespace property_tree {
template <typename Ch, typename Traits, typename E>
struct customize_stream<Ch, Traits, E*, void> {
  using stored_t = size_t;  // TODO convert to void* when it works
  static_assert(sizeof(stored_t) == sizeof(E*), "expected ptr width");

  static void insert(std::basic_ostream<Ch, Traits>& s, const E* e) {
    auto flags = s.flags();
    std::hex(s);  // write as hexadecimal
    std::showbase(s);
    s << reinterpret_cast<stored_t>(e);
    s.setf(flags);  // restore flags
  }
  static void extract(std::basic_istream<Ch, Traits>& s, E*& e) {
    auto flags = s.flags();
    std::hex(s);  // read as hexadecimal
    std::showbase(s);
    stored_t e_stored;
    s >> e_stored;
    s.setf(flags);  // restore flags
    e = reinterpret_cast<E*>(e_stored);
    if (!s.eof()) {
      s >> std::ws;
    }
  }
};
}  // namespace property_tree
}  // namespace boost

namespace mpqc {

namespace detail {
template <typename T>
struct register_keyval_ctor;

/// replaces boost serialization GUID: there is no default implementation of
/// guid()
/// to avoid the issues with the primary template getting picked up due to
/// incorrect
/// ordering of specialization and use
template <typename T>
const char* guid();

}  // namespace detail

class KeyVal;
namespace detail {
class SubTreeKeyVal;
}  // namespace detail

class DescribedClass;

template <typename T>
using Describable = std::is_base_of<DescribedClass, T>;

/// This class helps construction of (smart pointers to) C++ objects from ,
/// e.g., text input

/// To be able to construct class \c T from KeyVal objects do this:
/// <ol>
///   <li> make DescribedClass a public base of \c T, and </li>
///   <li> register \c T with DescribedClass . </li>
/// </ol>
/// The latter can be achieved in a number of ways, but the easiest is to add
/// any of the following statements
/// to a source file in the global scope:
/// <ol>
///   <li>if you want to use class name \c T as the type identifier in KeyVal
///   input: \c MPQC_CLASS_EXPORT_KEY(T) </li>
///   <li>if you want to use any other key \c Key as the type identifier in
///   KeyVal input: \c MPQC_CLASS_EXPORT_KEY2(T, Key) </li>
/// </ol>
/// It is the easiest to add these statements in the .cpp file that defines \c T
/// , but any other .cpp
/// file will work.
/// @note If \c T is a template class, you must register each instance of this
/// class you want to construct from KeyVal.
/// @warning to ensure that the class registration code of the derived class is
///          linked in, its destructor (at least) must be explicitly
///          instantiated.
///          Related: how gcc instantiates vtable and RTTI info see
///          <a
///          href="https://gcc.gnu.org/onlinedocs/gcc/Vague-Linkage.html">here</a>
/// @ingroup CoreKeyVal
class DescribedClass {
 public:
  DescribedClass() = default;
  virtual ~DescribedClass() {}

  ///
  typedef std::shared_ptr<DescribedClass> (*keyval_ctor_wrapper_type)(
      const KeyVal&);

  static keyval_ctor_wrapper_type type_to_keyval_ctor(
      const std::string& type_name);

  template <typename T>
  static void register_keyval_ctor() {
    const char* guid_T = mpqc::detail::guid<T>();
    assert(guid_T != nullptr);
    const std::string type_name(guid_T);
    auto& registry = keyval_ctor_registry_instance();
    assert(registry.find(type_name) == registry.end());
    registry.insert(std::make_pair(
        type_name,
        std::make_pair(keyval_ctor_wrapper<T>, std::cref(typeid(T)))));
  }

  /// query if T is registered
  /// @tparam T a class
  /// @return true if T is registered
  /// @warning cannot call this function before global object initialization has
  /// completed (i.e. after main() has started)
  template <typename T>
  static bool is_registered() {
    const auto& registry = keyval_ctor_registry_instance();
    for (const auto& elem : registry) {
      if (typeid(T) == elem.second.second) return true;
    }
    return false;
  }

  /// query the class key for T
  /// @tparam T a class
  /// @return class key with which T was registered, or empty string if T was
  /// not registered
  /// @warning cannot call this function before global object initialization has
  /// completed (i.e. after main() has started)
  template <typename T>
  static std::string class_key() {
    const auto& registry = keyval_ctor_registry_instance();
    for (const auto& elem : registry) {
      if (typeid(T) == elem.second.second) return elem.first;
    }
    return std::string();
  }

  /// query the class key for this object
  /// @return class key with which the class of this object was registered,
  ///         or empty string if it was not registered
  /// @warning cannot call this function before global object initialization has
  /// completed (i.e. after main() has started)
  std::string class_key() const {
    const auto& registry = keyval_ctor_registry_instance();
    for (const auto& elem : registry) {
      if (typeid(*this) == elem.second.second) return elem.first;
    }
    return std::string();
  }

  /// This class helps with registering DescribedClass with DescribedClass's
  /// static registry.

  /// To register the KeyVal ctor of type \c T create a single instance of this
  /// class
  /// @tparam T a class derived from DescribedClass
  /// @sa MPQC_CLASS_EXPORT_KEY2
  template <typename T, typename = std::enable_if_t<Describable<T>::value>>
  struct registrar {
    registrar() { DescribedClass::register_keyval_ctor<T>(); }
  };

 private:
  using keyval_ctor_registry_type =
      std::map<std::string,
               std::pair<keyval_ctor_wrapper_type,
                         std::reference_wrapper<const std::type_info>>>;

  // this is needed to force registry initialization BEFORE its use
  static keyval_ctor_registry_type& keyval_ctor_registry_instance() {
    static keyval_ctor_registry_type keyval_ctor_registry_;
    return keyval_ctor_registry_;
  }

  template <typename T>
  static std::shared_ptr<DescribedClass> keyval_ctor_wrapper(const KeyVal& kv) {
    return std::make_shared<T>(kv);
  }

 public:
  /// returns const ref to the keyval ctor registry
  static const keyval_ctor_registry_type& keyval_ctor_registry() {
    return const_cast<const keyval_ctor_registry_type&>(keyval_ctor_registry_instance());
  }

};
}  // namespace mpqc

/// @addtogroup CoreKeyVal
/// @{

/// MPQC_CLASS_REGISTER_GUID(K,T) associates key \c K as GUID for type \c T
/// \note this can be placed in a header or a source file.
#define MPQC_CLASS_REGISTER_GUID(K, ...) \
  namespace mpqc {                       \
  namespace detail {                     \
  template <>                            \
  const char* guid<__VA_ARGS__>() {      \
    return K;                            \
  }                                      \
  } /* detail */                         \
  } /* mpqc */                           \
/**/

/// MPQC_BOOST_CLASS_EXPORT_KEY2(K,T) associates key \c K with type \c T
/// \note this is a variadic version of BOOST_CLASS_EXPORT_KEY2 (see the docs
///       for Boost.Serialization); unlike BOOST_CLASS_EXPORT_KEY2 this can
///       register classes with two or more template arguments.
/// \note this should be placed in a header file (read
///       Boost.Serialization docs very carefully).
#define MPQC_BOOST_CLASS_EXPORT_KEY2(K, ...)               \
  namespace boost {                                        \
  namespace serialization {                                \
  template <>                                              \
  struct guid_defined<__VA_ARGS__> : boost::mpl::true_ {}; \
  namespace ext {                                          \
  template <>                                              \
  struct guid_impl<__VA_ARGS__> {                          \
    static inline const char* call() { return K; }         \
  };                                                       \
  }                                                        \
  } /* serialization */                                    \
  } /* boost */                                            \
/**/

/// \brief Associates a key (character string) with a class using
/// and register the class's KeyVal constructor with DescribedClass's registry.
/// This does not register class with Boost.Serialization, use
/// MPQC_BOOST_CLASS_EXPORT_KEY2 for that.
///
/// Use MPQC_CLASS_EXPORT_KEY2 to skip the KeyVal constructor
/// registration.
/// @note this should be placed in a source file, and only occur once for each
///       class.
#define MPQC_CLASS_EXPORT_KEY2(K, ...)                                         \
  MPQC_CLASS_REGISTER_GUID(K, __VA_ARGS__)                                     \
  namespace mpqc {                                                             \
  namespace detail {                                                           \
  template <>                                                                  \
  struct register_keyval_ctor<__VA_ARGS__> {                                   \
    static DescribedClass::registrar<__VA_ARGS__> const& r;                    \
  };                                                                           \
  DescribedClass::registrar<__VA_ARGS__> const&                                \
      register_keyval_ctor<__VA_ARGS__>::r =                                   \
          ::boost::serialization::singleton<                                   \
              DescribedClass::registrar<__VA_ARGS__>>::get_mutable_instance(); \
  }                                                                            \
  }                                                                            \
/**/

/// \brief Associates a key (character string) with a class
/// and register the class's KeyVal constructor with DescribedClass's registry.
/// This does not register class with Boost.Serialization, use
/// MPQC_BOOST_CLASS_EXPORT_KEY2 for that.
///
/// Identical to MPQC_CLASS_EXPORT_KEY2, but uses class name for the class key.
/// Use MPQC_BOOST_CLASS_EXPORT_KEY to skip the KeyVal ctor registration.
/// @sa MPQC_CLASS_EXPORT_KEY2
#define MPQC_CLASS_EXPORT_KEY(...) \
  MPQC_CLASS_EXPORT_KEY2(BOOST_PP_STRINGIZE(__VA_ARGS__), __VA_ARGS__)
/**/

/// \brief Forces the class instantiation so that it can be deserialized with
/// Boost.Serialization.
/// \note this is a variadic version of BOOST_CLASS_EXPORT_IMPLEMENT,
///       hence it can register classes with two or more template arguments.
/// \note this should be placed in a source file (read
///       Boost.Serialization docs very carefully).
#define MPQC_BOOST_CLASS_EXPORT_IMPLEMENT(...)                     \
  namespace boost {                                                \
  namespace archive {                                              \
  namespace detail {                                               \
  namespace extra_detail {                                         \
  template <>                                                      \
  struct init_guid<__VA_ARGS__> {                                  \
    static guid_initializer<__VA_ARGS__> const& g;                 \
  };                                                               \
  guid_initializer<__VA_ARGS__> const& init_guid<__VA_ARGS__>::g = \
      ::boost::serialization::singleton<                           \
          guid_initializer<__VA_ARGS__>>::get_mutable_instance()   \
          .export_guid();                                          \
  }                                                                \
  }                                                                \
  }                                                                \
  }                                                                \
  /**/

/// @}

namespace mpqc {

/**
    \brief KeyVal specifies C++ primitive data
    (booleans, integers, reals, string) and user-defined objects
    obtained from JSON/XML/INFO input or by programmatic construction.

    KeyVal is a (sub)tree of Key=Value pairs implemented with
    <a
   href="http://theboostcpplibraries.com/boost.propertytree">Boost.PropertyTree</a>.
    KeyVal extends the standard JSON/XML/INFO syntax to allow references as well
    as specification of registered C++ objects. See \ref keyval for the
   rationale,
    examples, and other details.

    KeyVal has reference semantics, i.e., copying KeyVal produces a KeyVal that
   refers to the same
    PropertyTree object and shares the class registry object with the original
   KeyVal object.

    @internal Since KeyVal is default-constructible and directly mutable, this
    obsoletes sc::AssignedKeyVal. Its behavior resembles PrefixKeyVal by
    combining prefix with a const tree. Hence this version of KeyVal roughly
    can be viewed as an assignable PrefixKeyVal of old MPQC, but supporting
    input from more modern formats like XML and JSON.

    @ingroup CoreKeyVal
 */
class KeyVal {
 public:
  /// data type for representing a property tree
  using ptree = boost::property_tree::iptree;
  using key_type = ptree::key_type;     // = std::string
  using value_type = ptree::data_type;  // = std::string
  constexpr static char separator = ':';

  /// reads the flag controlling whether reading a deprecated path results in an exception
  /// @return a boolean
  static bool throw_if_deprecated_path() {
    return throw_if_deprecated_path_accessor();
  }
  /// sets the flag controlling whether reading a deprecated path results in an exception
  /// @param f the new value of the flag
  static void set_throw_if_deprecated_path(bool f) {
    throw_if_deprecated_path_accessor() = f;
  }

 private:

  static bool& throw_if_deprecated_path_accessor() {
    static bool value = false;
    return value;
  }

  template <typename Class>
  struct is_sequence : std::false_type {};

  template <typename T, typename A>
  struct is_sequence<std::vector<T, A>> : std::true_type {};

  template <typename T, std::size_t N>
  struct is_sequence<std::array<T, N>> : std::true_type {};

  template <typename T, typename A>
  struct is_sequence<std::list<T, A>> : std::true_type {};

  template <typename T, std::size_t N>
  static void resize(std::array<T, N>& container, std::size_t size) {
    assert(size <= N);
  }

  template <typename T>
  static void resize(T& container, std::size_t size) {
    container.resize(size);
  }

  /// make_ptree is a helper to implement KeyVal::assign() for recursive data
  /// structures

  /// This specialization makes a ptree from a sequence.
  /// @tparam SequenceContainer any container for which
  /// KeyVal::is_sequence<SequenceContainer> is a \c std::true_type,
  ///         currently any of the following is allowed: \c std::array, \c
  ///         std::vector, \c std::list .
  /// @param value a sequence container to put at the path
  /// @param json_style if true, use empty keys so that JSON arrays are produced
  /// by KeyVal::write_json,
  ///                   else use 0-based integer keys, e.g. the first element
  ///                   will have key \c path:0,
  ///                   the second -- \c path:1, etc. (this is similar to how
  ///                   array elements were addressed in MPQC3)
  template <typename SequenceContainer>
  ptree make_ptree(
      const SequenceContainer& value, bool json_style = true,
      std::enable_if_t<KeyVal::is_sequence<SequenceContainer>::value>* =
          nullptr) {
    ptree result;
    size_t count = 0;
    for (const auto& v : value) {
      auto key =
          json_style
              ? key_type("")
              : std::to_string(count);  // assumes key_type == std::string
      ptree pt = make_ptree(v, json_style);
      result.push_back(std::make_pair(key, pt));
      ++count;
    }
    return result;
  }

  /// make_ptree is a helper to implement KeyVal::assign() for recursive data
  /// structures

  /// This specialization makes a ptree from an object of non-class non-sequence
  /// type.
  /// @tparam T a type which is not a sequence and not a class
  /// @param value a non-sequence object to put at the path
  template <typename T>
  std::enable_if_t<!KeyVal::is_sequence<T>::value &&
                       !utility::meta::can_construct<T, const KeyVal&>::value,
                   ptree>
  make_ptree(const T& value, bool) {
    ptree result;
    result.put("", value);
    return result;
  }

  /// counts the number of children of the node at this absolute path
  /// @param abs_path an absolute path
  /// @return an optional that does not contain a value if \c path does not exist or it points to a simple keyword,
  ///         otherwise containing the number of elements (>=0) if \c path points to an array or
  ///         a keyword group
  boost::optional<size_t> count_impl(const key_type& abs_path) const {
    auto child_opt =
        top_tree_->get_child_optional(ptree::path_type{abs_path, separator});
    if (child_opt == boost::optional<ptree&>())
      return boost::optional<size_t>{};
    else
      return child_opt->size();
  }

 public:
  /// creates empty KeyVal
  KeyVal();

  /// \brief copy ctor

  /// Since this class has reference semantics, the underlying data structures
  /// (ptree and class registry) are not copied.
  /// @note to clone top-level tree use KeyVal::clone()
  KeyVal(const KeyVal& other) = default;

  /// \brief copy assignment

  /// Since this class has reference semantics, the underlying data structures
  /// (ptree and class registry) are not copied.
  KeyVal& operator=(const KeyVal& other) = default;

  /// \brief creates a deep copy of this object
  /// \note the DescribedClass object registry is not copied
  KeyVal clone() const;

  /// \brief construct a KeyVal representing a subtree located at the given path

  /// Since this class has reference semantics, the result refers to the same
  /// underlying data structures (ptree and class registry).
  /// @internal corresponds to the old MPQC's PrefixKeyVal
  /// @param path the path to the subtree; absolute (within the top_tree) and
  ///        relative (with respect to this KeyVal's subtree) paths are allowed.
  /// @return the KeyVal object corresponding to \c path ;
  ///         if \c path does not exist, return an empty KeyVal
  virtual KeyVal keyval(const key_type& path) const {
    auto abs_path = resolve_path(path);
    return KeyVal(top_tree_, dc_registry_, abs_path);
  }

  /// \brief returns a shared_ptr to the (top) tree
  std::shared_ptr<ptree> top_tree() const { return top_tree_; }
  /// \brief returns a shared_ptr to this (sub)tree
  /// @note the result is aliased against the top tree
  std::shared_ptr<ptree> tree() const;

  /// checks whether the given path exists
  /// @param path the path
  /// @return true if \c path exists
  bool exists(const key_type& path) const {
    std::string resolved_path;
    bool bad_path = false;
    try {
      resolved_path = resolve_path(path);
    } catch (KeyVal::bad_input&) {
      bad_path = true;
    }
    return bad_path ? false : exists_(resolved_path);
  }

  /// check whether the given DescribedClass object exists in the registry already.
  /// @param path the path
  /// @return true if object \c path class exists
  /// TODO rename to exists_object_ptr
  bool exists_class(const key_type& path) const {
    bool exists_class_ptr = false;
    auto iter = dc_registry_->find(resolve_path(path));
    if (iter != dc_registry_->end()) {
      auto weak_ptr = iter->second;
      // if have unexpired cached ptr report true, otherwise false (and purge
      // the expired ptr)
      if (weak_ptr.expired())
        dc_registry_->erase(iter);
      else
        exists_class_ptr = true;
    }
    return exists_class_ptr;
  }

  /// counts the number of children of the node at this path
  /// @param path the path
  /// @return an optional that does not contain a value if \c path does not exist or it points to a simple keyword,
  ///         otherwise containing the number of elements (>=0) if \c path points to an array or
  ///         a keyword group
  boost::optional<size_t> count(const key_type& path) const {
    auto abs_path = resolve_path(path);
    return count_impl(abs_path);
  }

  /// @brief assign simple \c value at the given path (overwrite, if necessary)
  /// @param value a "simple" value, i.e. it can be converted to
  /// a KeyVal::key_type using a
  /// std::basic_ostream<KeyVal::key_type::value_type>
  template <typename T,
            typename = std::enable_if_t<not KeyVal::is_sequence<T>::value>>
  KeyVal& assign(const key_type& path, const T& value) {
    auto abs_path = to_absolute_path(path);
    top_tree_->put(ptree::path_type{abs_path, separator}, value);
    return *this;
  }

  /// @brief Assign a sequence container at the given path (overwrite, if
  /// necessary).

  /// @tparam SequenceContainer any container for which
  /// KeyVal::is_sequence<SequenceContainer> is a \c std::true_type,
  ///         currently any of the following is allowed: \c std::array, \c
  ///         std::vector, \c std::list . Can be a recursive sequence, i.e.
  ///         sequence of sequence.
  /// @param path the target path
  /// @param value a sequence container to put at the path
  /// @param json_style if true, use empty keys so that JSON arrays are produced
  /// by KeyVal::write_json,
  ///                   else use 0-based integer keys, e.g. the first element
  ///                   will have key \c path:0,
  ///                   the second -- \c path:1, etc. (this is similar to how
  ///                   array elements were addressed in MPQC3)
  template <typename SequenceContainer>
  KeyVal& assign(
      const key_type& path, const SequenceContainer& value,
      bool json_style = true,
      std::enable_if_t<KeyVal::is_sequence<SequenceContainer>::value>* =
          nullptr) {
    auto abs_path = to_absolute_path(path);
    ptree obj = make_ptree(value, json_style);
    top_tree_->add_child(ptree::path_type{abs_path, separator}, obj);
    return *this;
  }

  /// assign the given pointer to a DescribedClass object at the given path
  /// (overwrite, if necessary)
  /// @tparam T a class directly derived from DescribedClass
  /// @param path the path to \c value
  /// @param value the object pointer to assign to path \c path
  /// @warning these key/value pairs are not part of ptree, hence cannot be
  /// written to JSON/XML
  template <typename T = DescribedClass,
            typename = std::enable_if_t<Describable<T>::value>>
  KeyVal& assign(const key_type& path, const std::shared_ptr<T>& value) {
    auto dc_value = std::static_pointer_cast<DescribedClass>(value);
    std::weak_ptr<DescribedClass> weak_value = dc_value;
    auto abs_path = to_absolute_path(path);
    (*dc_registry_)[abs_path] = weak_value;
    return *this;
  }

  /// erases entries located under \c path
  /// @param path the target path
  KeyVal& erase(const key_type& path) {
    auto abs_path = resolve_path(path);
    key_type abs_path_basename, name;
    std::tie(abs_path_basename, name) = path_basename(abs_path);
    top_tree_->get_child(ptree::path_type{abs_path_basename, separator})
        .erase(name);  // count is ignored
    return *this;
  }

  // validators
  struct dummy_validator_t {
    template <typename... Args>
    constexpr bool operator()(Args&&... args) const {
      return true;
    }
  };
  static constexpr dummy_validator_t dummy_validator{};
  struct is_nonnegative_t {
    template <typename Arg>
    constexpr bool operator()(Arg&& arg) const {
      return arg >= 0;
    }
  };
  static const is_nonnegative_t is_nonnegative;
  struct is_positive_t {
    template <typename Arg>
    constexpr bool operator()(Arg&& arg) const {
      return arg > 0;
    }
  };
  static const is_positive_t is_positive;

  /// return the result of converting the value at the given path to the desired
  /// type.

  /// @tparam T the desired value type, must be a "simple" type that can be
  /// accepted by KeyVal::assign()
  /// @tparam Validator the type of @c validator
  /// @param path the path to the value
  /// @param validator a callable for which @c validator(T) returns a boolean.
  ///        Before returning a value, it will be validated by @c validator .
  ///        The default is a dummy validator that always returns @c true.
  /// @return value stored at \c path converted to type \c T
  /// @throws KeyVal::bad_input if path not found or cannot convert value
  /// representation to the desired type
  /// @throws KeyVal::bad_input if validation failed.
  template <typename T, typename Validator = dummy_validator_t,
            typename = std::enable_if_t<
                not KeyVal::is_sequence<T>::value &&
                not utility::meta::can_construct<T, const KeyVal&>::value &&
                std::is_same<bool, decltype(std::declval<Validator>()(
                                       std::declval<const T&>()))>::value>>
  T value(const key_type& path, Validator&& validator = Validator{}) const {
    T result;
    try {
      if (auto subtree = this->get_subtree(path)) {
        result = subtree.get().get_value<T>();
      } else
        throw KeyVal::bad_input(std::string("path not found"),
                                to_absolute_path(path));
    } catch (boost::property_tree::ptree_bad_data& x) {
      throw KeyVal::bad_input(
          std::string("value ") + x.data<KeyVal::value_type>() +
              " could not be converted to the desired datatype",
          to_absolute_path(path));
    }
    if (validator(result) == false) {
      std::ostringstream oss;
      oss << "value " << result << " did not pass validation";
      throw KeyVal::bad_input(oss.str(), to_absolute_path(path));
    }
    return result;
  }

  /// return the result of converting the value at the given path to the desired
  /// type.

  /// @tparam T the desired value type, must be a "simple" type that can be
  /// accepted by KeyVal::assign()
  /// @tparam Validator the type of @c validator
  /// @param path the path to the value
  /// @param validator a callable for which @c validator(T) returns a boolean.
  ///        Before returning a value, it will be validated by @c validator .
  ///        The default is a dummy validator that always returns @c true.
  /// @return value stored at \c path converted to type \c T
  /// @throws KeyVal::bad_input if path not found or cannot convert value
  /// @throws KeyVal::bad_input if validation failed.
  /// representation to the desired type
  template <typename T, typename Validator = dummy_validator_t>
  std::enable_if_t<not KeyVal::is_sequence<T>::value &&
                       utility::meta::can_construct<T, const KeyVal&>::value &&
                       std::is_same<bool, decltype(std::declval<Validator>()(
                                              std::declval<const T&>()))>::value,
                   T>
  value(const key_type& path, Validator&& validator = Validator{}) const;

  /// return value corresponding to a path and convert to a std::vector.
  /// @tparam T the desired value type
  /// @tparam Validator the type of @c validator
  /// @param path the path
  /// @param validator a callable for which @c validator(T) returns a boolean.
  ///        Before returning a value, it will be validated by @c validator .
  ///        The default is a dummy validator that always returns @c true.
  /// @throw KeyVal::bad_input if \c path is bad
  /// @throws KeyVal::bad_input if validation failed.
  /// @return value of type \c T
  template <typename SequenceContainer, typename Validator = dummy_validator_t>
  SequenceContainer value(
      const key_type& path, Validator&& validator = Validator{},
      std::enable_if_t<
          KeyVal::is_sequence<SequenceContainer>::value &&
          std::is_same<bool, decltype(std::declval<Validator>()(
                                 std::declval<const typename SequenceContainer::
                                                  value_type&>()))>::value>* =
          nullptr) const;  // implemented after SubTreeKeyval is implemented

  /// return value corresponding to a path and convert to the desired type.
  /// @tparam T the desired value type
  /// @tparam Validator the type of @c validator
  /// @param path the path to the value
  /// @param default_value the default value will be returned if @c path is not
  /// found.
  /// @param deprecated_path a deprecated path will only be queried if its not
  ///        empty and @c path not found;
  ///        if its value is used a message will be added
  ///        to @std::cerr. The default is an empty string.
  /// @param validator a callable for which @c validator(T) returns a boolean.
  ///        Before returning a value, it will be validated by @c validator .
  ///        The default is a dummy validator that always returns @c true.
  /// @return value of type \c T
  /// @throws KeyVal::bad_input if validation failed.
  template <typename T, typename Validator = dummy_validator_t,
            typename = std::enable_if_t<
                not KeyVal::is_sequence<T>::value &&
                std::is_same<bool, decltype(std::declval<Validator>()(
                                       std::declval<const T&>()))>::value>>
  T value(const key_type& path, const T& default_value,
          const key_type& deprecated_path,
          Validator&& validator = Validator{}) const {
    T result = default_value;
    const key_type* read_path = &path;

    if (auto subtree = this->get_subtree(path)) {
      result = subtree.get().template get_value<T>(default_value);
    } else if (!deprecated_path.empty() &&
               (subtree = this->get_subtree(deprecated_path))) {
      read_path = &deprecated_path;
      auto result_optional = subtree.get().template get_value_optional<T>();
      if (result_optional) {
        result = result_optional.get();
        if (KeyVal::throw_if_deprecated_path())
          throw ;
        else
          std::cerr << "KeyVal read value from deprecated path "
                    << deprecated_path << std::endl;
      }
    }

    if (validator(result) == false) {
      std::ostringstream oss;
      oss << "value " << result << " did not pass validation";
      throw KeyVal::bad_input(oss.str(), to_absolute_path(*read_path));
    }
    return result;
  }

  /// return value corresponding to a path and convert to the desired type.
  /// @tparam T the desired value type
  /// @tparam Validator the type of @c validator
  /// @param path the path to the value
  /// @param default_value the default value will be returned if @c path is not
  /// found.
  /// @param validator a callable for which @c validator(T) returns a boolean.
  ///        Before returning a value, it will be validated by @c validator .
  ///        The default is a dummy validator that always returns @c true.
  /// @return value of type \c T
  /// @throws KeyVal::bad_input if validation failed.
  template <typename T, typename Validator = dummy_validator_t,
            typename = std::enable_if_t<
                not KeyVal::is_sequence<T>::value &&
                std::is_same<bool, decltype(std::declval<Validator>()(
                                       std::declval<const T&>()))>::value>>
  T value(const key_type& path, const T& default_value,
          Validator&& validator = Validator{}) const {
    T result = default_value;

    if (auto subtree = this->get_subtree(path)) {
      result = subtree.get().template get_value<T>(default_value);
    }

    if (validator(result) == false) {
      std::ostringstream oss;
      oss << "value " << result << " did not pass validation";
      throw KeyVal::bad_input(oss.str(), to_absolute_path(path));
    }
    return result;
  }

  /// return a pointer to an object at the given path

  /// Constructs a smart pointer to object of type \c T . If user provided
  /// keyword "type" in this KeyVal object its value will be used for the object
  /// type.
  /// Otherwise, if \c T is a <i>concrete</i> type,
  /// this will construct object of type \c T. In either case, the type of
  /// constructed object must have DescribedClass as its public base <i>and</i>
  /// registered.
  /// @tparam T a class derived from DescribedClass
  /// @param path the path; if not provided, use the empty path
  /// @param bypass_registry if \c true will not query the registry for the
  /// existence of the object
  ///        and will not place the newly-constructed object in the registry
  /// @return a std::shared_ptr to \c T ; null pointer will be returned if \c
  /// path does not exist.
  /// @throws KeyVal::bad_input if the value of keyword "type" is not a
  /// registered class key,
  /// or if the user did not specify keyword "type" and class T is abstract or
  /// is not registered.
  /// @note `class_ptr<T>("path")`  is equivalent to
  /// `keyval("path").class_ptr<T>()`
  /// TODO rename to object_ptr()
  template <typename T = DescribedClass,
            typename = std::enable_if_t<Describable<T>::value>>
  std::shared_ptr<T> class_ptr(const key_type& path = key_type(),
                               bool bypass_registry = false) const {
    // if this class already exists in the registry under path
    // (e.g. if the ptr was assigned programmatically), return immediately
    if (!bypass_registry) {
      auto it = dc_registry_->find(path);
      if (it != dc_registry_->end()) {
        // if have unexpired weak ptr, convert to shared and return
        auto weak_ptr = it->second;
        if (!weak_ptr.expired())
          return std::dynamic_pointer_cast<T>(weak_ptr.lock());
      }
    }

    // otherwise, resolve the path and build the class (or pick up the cached
    // copy)
    auto abs_path = resolve_path(path);

    // if this class already exists, return the ptr
    if (!bypass_registry) {
      auto it = dc_registry_->find(abs_path);
      if (it != dc_registry_->end()) {
        // if have unexpired weak ptr, convert to shared and return
        auto weak_ptr = it->second;
        if (!weak_ptr.expired())
          return std::dynamic_pointer_cast<T>(weak_ptr.lock());
      }
    }

    // return nullptr if the path does not exist, or does not point to an aggregate
    if (!this->exists_(abs_path) || !static_cast<bool>(this->count_impl(abs_path))) {
      return std::shared_ptr<T>();
    }

    // compute the type name of the constructed object
    std::string result_type_name;
    auto type_path = concat_path(abs_path, "type");
    if (not exists_(type_path)) {
      // no user override for type = try to construct T
      if (std::is_abstract<T>::value) {
        throw KeyVal::bad_input(
            std::string(
                "KeyVal::class_ptr<T>(): T is an abstract class, "
                "hence keyword \"type\" must be given for object at path " +
                path),
            abs_path);
      }
      if (!DescribedClass::is_registered<T>()) {
        throw KeyVal::bad_input(
            std::string(
                "KeyVal::class_ptr<T>(): T is an unregistered concrete class, "
                "hence keyword \"type\" must be given for object at path " +
                path),
            abs_path);
      }
      result_type_name = DescribedClass::class_key<T>();
    } else {  // user provided \"type\" keyword, use it to override the type
      result_type_name =
          top_tree_->get<std::string>(ptree::path_type{type_path, separator});
    }

    // map class type to its ctor
    auto ctor = DescribedClass::type_to_keyval_ctor(result_type_name);
    // call the ctor
    // if nonempty path, construct a KeyVal object rooted at that path,
    // otherwise use self
    auto result = !path.empty() ? (*ctor)(this->keyval(path)) : (*ctor)(*this);
    if (!bypass_registry) (*dc_registry_)[abs_path] = result;
    return std::dynamic_pointer_cast<T>(result);
  }

  /// read/write KeyVal specified in one of accepted Boost.PropertyTree formats
  /// (see
  /// <a
  /// href="http://www.boost.org/doc/libs/master/doc/html/property_tree/parsers.html">Boost.PropertyTree
  /// docs</a>)
  /// @{

  /// write to stream in XML form
  void write_xml(std::basic_ostream<typename key_type::value_type>& os) const {
    boost::property_tree::xml_parser::write_xml(os, *(tree()));
  }
  /// write from stream in XML form
  void read_xml(std::basic_istream<typename key_type::value_type>& is) {
    boost::property_tree::xml_parser::read_xml(is, *(tree()));
  }

  /// write to stream in JSON form
  void write_json(std::basic_ostream<key_type::value_type>& os) const {
    boost::property_tree::json_parser::write_json(os, *(tree()));
  }
  /// read from stream in JSON form
  void read_json(std::basic_istream<key_type::value_type>& is) {
    boost::property_tree::json_parser::read_json(is, *(tree()));
  }

  /// write to stream in INFO form
  void write_info(std::basic_ostream<typename key_type::value_type>& os) const {
    boost::property_tree::info_parser::write_info(os, *(tree()));
  }
  /// write from stream in INFO form
  void read_info(std::basic_istream<typename key_type::value_type>& is) {
    boost::property_tree::info_parser::read_info(is, *(tree()));
  }
  /// @}

  /// @return the path from the root of top tree to this subtree
  std::string path() const { return path_; }

 private:
  std::shared_ptr<ptree> top_tree_;
  // 'dc' = DescribedClass
  using dc_registry_type = std::map<std::string, std::weak_ptr<DescribedClass>>;
  std::shared_ptr<dc_registry_type> dc_registry_;
  const key_type path_;  //!< path from the top of \c top_tree_ to this subtree

  /// creates a KeyVal
  KeyVal(const std::shared_ptr<ptree>& top_tree,
         const std::shared_ptr<dc_registry_type>& dc_registry,
         const key_type& path)
      : top_tree_(top_tree), dc_registry_(dc_registry), path_(path) {}

  /// creates a KeyVal using am existing \c ptree
  /// @param pt a ptree
  KeyVal(ptree pt)
      : top_tree_(std::make_shared<ptree>(std::move(pt))),
        dc_registry_(),
        path_("") {}

  /// given a path that contains \c ".." elements, returns the equivalent path
  /// without such elements. For example: \c realpath("tmp:..:x") returns \c "x"
  /// @param path the path optionally including \c ".."
  /// @return the path with the \c ".." elements resolved
  /// @note this does not resolve references or convert to absolute path
  /// @throw KeyVal::bad_input if path is invalid; an example is ".." .
  static key_type realpath(const key_type& path) {
    key_type result = path;
    const key_type parent_token("..");
    auto parent_token_location = result.find(parent_token);
    while (parent_token_location != key_type::npos) {
      // make sure .. is preceeded by :
      if (parent_token_location == 0)
        throw KeyVal::bad_input(
            "KeyVal: reference path cannot begin with \"..\"", path);
      if (path[parent_token_location - 1] != separator)
        throw KeyVal::bad_input(
            std::string(
                "KeyVal: \"..\" in reference path must be preceeded by \"") +
                separator + "\"",
            path);
      // did we find '..:' or '..\0'?
      auto parent_token_length =
          (result.find(separator, parent_token_location + 1) ==
           parent_token_location + 2)
              ? 3
              : 2;
      // find the beginning of the preceding path element
      // NB handle the corner case of path starting with ":.."
      auto preceding_token_location =
          parent_token_location > 1
              ? path.rfind(separator, parent_token_location - 2)
              : 0;
      if (preceding_token_location == key_type::npos)  // if separator not
                                                       // found, the preceeding
                                                       // element is the
                                                       // beginning of the path
        preceding_token_location = 0;
      result.erase(preceding_token_location,
                   (parent_token_location + parent_token_length -
                    preceding_token_location));
      parent_token_location = result.find(parent_token);
    }
    return result;
  }

  /// computes an absolute path from a given path prefix and an (absolute or
  /// relative) path
  /// @param path_prefix
  /// @param path an absolute path (starts with "$:") or a relative path (the
  /// leading "$" is dropped)
  /// @return the absolute path
  /// @note this does not resolve references
  /// @throw KeyVal::bad_input if path is invalid; an example is ".." .
  static key_type to_absolute_path(const key_type& path_prefix,
                                   const key_type& path) {
    auto is_ref = path.size() != 0 && path[0] == '$';
    auto is_abs_ref = is_ref && path.size() > 1 && path[1] == ':';
    key_type abs_path;
    key_type path_copy = path;
    if (is_abs_ref) {
      abs_path = path_copy.substr(2);  // drop the leading "$:"
    } else {
      key_type rel_path =
          is_ref ? path_copy.substr(1) : path_copy;  // drop the leading "$"
      abs_path = concat_path(path_prefix, rel_path);
    }
    abs_path = realpath(abs_path);
    return abs_path;
  }

  /// computes an absolute path from a given (absolute or relative) path
  /// @param path an absolute path (starts with "$:") or a relative path (the
  /// leading "$" is dropped)
  /// @return the absolute path
  /// @note this does not resolve references
  /// @throw KeyVal::bad_input if path is invalid; an example is ".." .
  key_type to_absolute_path(const key_type& path) const {
    return to_absolute_path(path_, path);
  }

  /// Concatenates two paths
  /// @return a path obtained by concatenating \c prefix and \c postfix ;
  static key_type concat_path(const key_type& prefix, const key_type& postfix) {
    // if (postfix.size() > 0) assert(postfix[0] != separator);
    auto result =
        (postfix == "")
            ? prefix
            : ((prefix == "") ? postfix : prefix + separator + postfix);
    return result;
  }

  /// normalizes path by 1. converting path to absolute path (see \c
  /// to_absolute_path())
  /// and 2. resolving any refs in the path
  /// @throw KeyVal::bad_input if path is invalid; an example is ".." .
  key_type resolve_path(const key_type& path) const {
    auto abs_path = to_absolute_path(path);
    auto result = resolve_refs(abs_path);
    return result;
  }

  /// This resolves the references in the path.
  /// @note repeatedly recurses, hence can be expensive
  /// @note circular references are not detected, thus to prevent infinite
  /// recursion will throw if too many iterations are needed
  /// @param niter_max max number of iterations
  /// @throws \c KeyVal::bad_input if max number of iterations exceeded
  key_type resolve_refs(const key_type& path, size_t niter_max = 10) const {
    key_type result = path;
    size_t niter = 0;
    while (niter < niter_max) {
      bool done = resolve_first_ref(result);
      if (done) return result;
      ++niter;
    }
    // too many iterations if still here
    throw KeyVal::bad_input("excessive or circular references in path", path);
  }

  /// @returns true if \c path did not include a reference, even if path does
  /// not exist
  bool resolve_first_ref(key_type& path) const {
    if (path.size() == 0)
      return true;  // handle this corner case now to make the rest a bit
                    // cleaner

    auto subpath_last_char = path.find(
        separator, 1);  // if path starts ":" no need to check root of the tree
    if (subpath_last_char == key_type::npos)
      subpath_last_char =
          path.size();  // this ensures that the full path is checked
    while (subpath_last_char <= path.size()) {
      auto subpath = path.substr(0, subpath_last_char);
      if (!exists_(subpath)) return true;
      auto value = top_tree_->get<key_type>(ptree::path_type{
          subpath, separator});  // value corresponding to subpath
      auto subpath_is_a_ref = (value.size() > 0 && value[0] == '$');
      if (subpath_is_a_ref) {
        // resolve subpath key + value to a new subpath
        auto subpath_base = path_pop(subpath);
        auto new_subpath = to_absolute_path(subpath_base, value);

        // replace the old subpath with the new subpath
        path = new_subpath + path.substr(subpath.size());

        return false;
      }
      if (subpath_last_char ==
          path.size())  // subpath = path and it's not a ref? done
        return true;
      subpath_last_char = path.find(separator, subpath_last_char + 1);
      if (subpath_last_char == key_type::npos)
        subpath_last_char =
            path.size();  // this ensures that the full path is checked
    }
    return true;  // no refs found
  }

  /// pops off last ":token" (if separator not found, the result is an empty
  /// path)
  static key_type path_pop(const key_type& path) {
    key_type result = path;
    auto last_separator_location = result.rfind(separator);
    if (last_separator_location == key_type::npos) last_separator_location = 0;
    result.erase(last_separator_location);
    return result;
  }

  /// returns {basename,key} for the given path; e.g.
  /// for path ":key1:key2:key3:key4" returns {":key1:key2:key3", "key4"}
  static std::tuple<key_type, key_type> path_basename(const key_type& path) {
    auto last_separator_location = path.rfind(separator);
    if (last_separator_location == key_type::npos) last_separator_location = 0;
    key_type basename = path;
    basename.erase(last_separator_location);
    key_type key = (last_separator_location == 0)
                       ? path
                       : path.substr(last_separator_location + 1);
    return std::make_tuple(basename, key);
  }

  /// checks whether the given path exists; does not resolve path unlike
  /// KeyVal::exists()
  /// @param path the path
  /// @return true if \c path exists
  bool exists_(const key_type& path) const {
    return top_tree_->get_child_optional(ptree::path_type{path, separator}) !=
           boost::optional<ptree&>();
  }

 protected:
  /// returns a reference to a subtree located at (absolute or relative) path.
  /// @param path the path to the subtree
  virtual boost::optional<const ptree&> get_subtree(
      const key_type& path) const {
    auto abs_path = resolve_path(path);
    return const_cast<const ptree&>(*top_tree_)
        .get_child_optional(ptree::path_type{abs_path, separator});
  }

 public:
  /// \brief KeyVal exception
  struct bad_input : public std::runtime_error {
    bad_input(const std::string& _what, const key_type& path)
        : std::runtime_error(_what + " (path=" + path + ")") {}
    virtual ~bad_input() noexcept {}
  };
};  // KeyVal

/// union of two KeyVal objects
/// @note obsoletes sc::AggregateKeyVal
/// @ingroup KeyValCore
KeyVal operator+(const KeyVal& first, const KeyVal& second);

namespace detail {
/// is a KeyVal that instead of keeping a path to the root of this subtree has
/// actual ptree ref
/// currently only useful for parsing JSON arrays (with their braindead empty
/// element keys),
/// specifically is mandatory for arrays of objects
class SubTreeKeyVal : public KeyVal {
 public:
  SubTreeKeyVal(const ptree& tree, const KeyVal& kv)
      : KeyVal(kv), subtree_(tree) {}

  KeyVal keyval(const key_type& path) const override {
    if (auto subtree = this->get_subtree(path))
      return SubTreeKeyVal(subtree.get(), *this);
    else
      throw KeyVal::bad_input(
          std::string("detail::KeyVal::keyval: path not found"), path);
  }

 private:
  const ptree& subtree_;

  /// overrides KeyVal::
  boost::optional<const ptree&> get_subtree(
      const key_type& path) const override {
    if (path.empty())  // getting elements of arrays of scalars with give empty
                       // path
      return subtree_;
    else if (path[0] != separator)  // will parse down into the subtree when
                                    // making arrays of objects
      return subtree_.get_child_optional(ptree::path_type{path, separator});
    else  // what should the default be if given absolute path? option 1: use
          // KeyVal::get_subtree,
      // else throw, probably did not want SubTreeKeyVal in this case after all
      throw KeyVal::bad_input("detail::SubTreeKeyVal given an absolute path",
                              path);
    // return this->KeyVal::get_subtree(path);
  }
};
}  // namespace detail

template <typename SequenceContainer, typename Validator>
SequenceContainer KeyVal::value(
    const key_type& path,
    Validator&& validator_,
    std::enable_if_t<KeyVal::is_sequence<SequenceContainer>::value &&
        std::is_same<bool, decltype(std::declval<Validator>()(
            std::declval<const typename SequenceContainer::
            value_type&>()))>::value>*) const {
  using value_type = typename SequenceContainer::value_type;
  SequenceContainer result;
  std::remove_reference_t<Validator> validator(std::forward<Validator>(validator_));
  if (auto vec_ptree_opt = this->get_subtree(path)) {
    try {
      auto vec_ptree = vec_ptree_opt.get();
      KeyVal::resize(result, vec_ptree.size());
      auto iter = result.begin();
      size_t count = 0;
      for (const auto& elem_ptree : vec_ptree) {
        assert(elem_ptree.first == ""  // JSON array spec
               ||
               elem_ptree.first ==
                   std::to_string(
                       count)  // 0-based array keys, a la ipv2, usable with XML
        );
        detail::SubTreeKeyVal stree_kv(elem_ptree.second, *this);
        *iter = stree_kv.value<value_type>("", validator);
        ++iter;
        ++count;
      }
      return result;
    } catch (boost::property_tree::ptree_bad_data& x) {
      throw KeyVal::bad_input(
          std::string("value ") + x.data<KeyVal::value_type>() +
              " could not be converted to the desired datatype",
          resolve_path(path));
    } catch (boost::property_tree::ptree_bad_path& x) {
      throw KeyVal::bad_input(
          std::string("path ") + x.path<key_type>() + " not found",
          resolve_path(path));
    }
  } else
    throw KeyVal::bad_input(std::string("path not found"),
                            to_absolute_path(path));
}

template <typename T, typename Validator>
std::enable_if_t<not KeyVal::is_sequence<T>::value &&
                     utility::meta::can_construct<T, const KeyVal&>::value &&
    std::is_same<bool, decltype(std::declval<Validator>()(
        std::declval<const T&>()))>::value,
                 T>
KeyVal::value(const key_type& path, Validator&& validator) const {
  const detail::SubTreeKeyVal stree_kv(this->get_subtree(path).get(), *this);
  auto result = T(stree_kv);
  if (validator(result) == false) {
    std::ostringstream oss;
    oss << "value " << utility::to_ostream(result) << " did not pass validation";
    throw KeyVal::bad_input(oss.str(), to_absolute_path(path));
  }
  return result;
}

static_assert(utility::meta::can_construct<double, const KeyVal&>::value ==
                  false,
              "");
static_assert(utility::meta::can_construct<std::string, const KeyVal&>::value ==
                  false,
              "");

}  // namespace mpqc

/// @}

#endif  // MPQC4_SRC_MPQC_UTIL_KEYVAL_KEYVAL_H_
