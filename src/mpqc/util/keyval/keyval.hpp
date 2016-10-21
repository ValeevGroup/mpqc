/*
 * keyval.hpp
 *
 *  Created on: Apr 18, 2016
 *      Author: evaleev
 */

#ifndef SRC_MPQC_UTIL_KEYVAL_KEYVAL_HPP_
#define SRC_MPQC_UTIL_KEYVAL_KEYVAL_HPP_

#include <cassert>
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <array>
#include <list>
#include <stdexcept>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/serialization/export.hpp>

// serialize all pointers as void*
// NB XCode 7.3.1 (7D1014) libc++ char stream does not properly deserialize
// void*, use size_t instead
namespace boost { namespace property_tree
{
    template <typename Ch, typename Traits, typename E>
    struct customize_stream<Ch,Traits,E*,void>
    {
        using stored_t  = size_t; // TODO convert to void* when it works
        static_assert(sizeof(stored_t) == sizeof(E*), "expected ptr width");

        static void insert(std::basic_ostream<Ch, Traits>& s, const E* e) {
            auto flags = s.flags();
            std::hex(s); // write as hexadecimal
            std::showbase(s);
            s << reinterpret_cast<stored_t>(e);
            s.setf(flags); // restore flags
        }
        static void extract(std::basic_istream<Ch, Traits>& s, E*& e) {
            auto flags = s.flags();
            std::hex(s); // read as hexadecimal
            std::showbase(s);
            stored_t e_stored;
            s >> e_stored;
            s.setf(flags); // restore flags
            e = reinterpret_cast<E*>(e_stored);
            if(!s.eof()) {
                s >> std::ws;
            }
        }
    };
}} // namespace boost::property_tree


#include "mpqc/util/misc/type_traits.h"

namespace mpqc {

class KeyVal;

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
///          linked in its destructor (at least) must be explicitly instantiated
/// @ingroup CoreKeyVal
class DescribedClass {
 public:
  DescribedClass() = default;
  virtual ~DescribedClass() = default;

  ///
  typedef std::shared_ptr<DescribedClass> (*keyval_ctor_wrapper_type)(
      const KeyVal&);

  static keyval_ctor_wrapper_type type_to_keyval_ctor(
      const std::string& type_name);

  template <typename T>
  static void register_keyval_ctor() {
    const std::string type_name = boost::serialization::guid<T>();
    auto& registry = keyval_ctor_registry();
    assert(registry.find(type_name) == registry.end());
    registry[type_name] = keyval_ctor_wrapper<T>;
  }

  /// This class helps with registering DescribedClass with DescribedClass's
  /// static registry.

  /// To register the KeyVal ctor of type \c T create a single instance of this
  /// class
  /// @tparam T a class derived from DescribedClass
  /// @sa MPQC_CLASS_EXPORT_KEY2
  template <typename T, typename = enable_if_t<Describable<T>::value>>
  struct registrar {
    registrar() { DescribedClass::register_keyval_ctor<T>(); }
  };

 private:
  using keyval_ctor_registry_type =
      std::map<std::string, keyval_ctor_wrapper_type>;

  // this is needed to force registry initialization BEFORE its use
  static keyval_ctor_registry_type& keyval_ctor_registry() {
    static keyval_ctor_registry_type keyval_ctor_registry_;
    return keyval_ctor_registry_;
  }

  template <typename T>
  static std::shared_ptr<DescribedClass> keyval_ctor_wrapper(const KeyVal& kv) {
    return std::make_shared<T>(kv);
  }
};
}  // namespace mpqc

namespace mpqc {
namespace detail {
template <typename T>
struct register_keyval_ctor;
}
}

/// @addtogroup CoreKeyVal
/// @{

/// MPQC_BOOST_CLASS_EXPORT_KEY2(K,T) associates key \c K with type \c T
/// \note this is a variadic version of BOOST_CLASS_EXPORT_KEY2
#define MPQC_BOOST_CLASS_EXPORT_KEY2(K, ...)               \
  namespace boost {                                        \
  namespace serialization {                                \
  template <>                                              \
  struct guid_defined<__VA_ARGS__> : boost::mpl::true_ {}; \
  template <>                                              \
  inline const char* guid<__VA_ARGS__>() {                 \
    return K;                                              \
  }                                                        \
  } /* serialization */                                    \
  } /* boost */                                            \
/**/

/// \brief Associates a key (character string) with a class using
/// Boost.Serialization
/// and register the class's KeyVal constructor with DescribedClass's registry.
///
/// Use MPQC_BOOST_CLASS_EXPORT_KEY2 to skip the KeyVal constructor
/// registration.
#define MPQC_CLASS_EXPORT_KEY2(K, ...)                                         \
  MPQC_BOOST_CLASS_EXPORT_KEY2(K, __VA_ARGS__)                                 \
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

/// \brief Associates a key (character string) with a class using
/// Boost.Serialization
/// and register the class's KeyVal constructor with DescribedClass's registry.
///
/// Identical to MPQC_CLASS_EXPORT_KEY2, but uses class name for the class key.
/// Use MPQC_BOOST_CLASS_EXPORT_KEY to skip the KeyVal ctor registration.
/// @sa MPQC_CLASS_EXPORT_KEY2
#define MPQC_CLASS_EXPORT_KEY(...) \
  MPQC_CLASS_EXPORT_KEY2(BOOST_PP_STRINGIZE(__VA__ARGS__), __VA_ARGS__)
/**/

/// \brief Forces the class instantiation so that it can be deserialized with
/// Boost.Serialization and/or constructed from a KeyVal.
/// \note this is a variadic version of BOOST_CLASS_EXPORT_IMPLEMENT
#define MPQC_CLASS_EXPORT_IMPLEMENT(...)                           \
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
    obtained from JSON/XML input or by programmatic construction.

    KeyVal is a (sub)tree of Key=Value pairs implemented with
    <a
   href="http://theboostcpplibraries.com/boost.propertytree">Boost.PropertyTree</a>.
    KeyVal extends the standard JSON/XML syntax to allow references as well as
    specification of registered C++ objects. See \ref keyval for the rationale,
   examples,
    and other details.

    @note Since KeyVal is default-constructible and directly mutable, this
    obsoletes sc::AssignedKeyVal. Its behavior resembles PrefixKeyVal by
    combining prefix with a const tree. Hence this version of KeyVal roughly
    can be viewed as an assignable PrefixKeyVal of old MPQC, but supporting
    input from more modern formats like XML and JSON.

    @ingroup CoreKeyVal
 */
class KeyVal {
 private:
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

 public:
  /// data type for representing a property tree
  using ptree = boost::property_tree::iptree;
  using key_type = ptree::key_type;     // = std::string
  using value_type = ptree::data_type;  // = std::string
  constexpr static char separator = ':';

  /// creates empty KeyVal
  KeyVal();

  /// (mostly shallow) copy ctor (to clone top-level tree use top_tree())
  KeyVal(const KeyVal& other) = default;
  /// copy assignment
  KeyVal& operator=(const KeyVal& other) = default;

  /// construct a KeyVal representing a subtree located at the given path
  /// @note corresponds to the old MPQC's PrefixKeyVal
  /// @param path the path to the subtree; absolute (within the top_tree) and
  ///        relative (with respect to this KeyVal's subtree) paths are allowed.
  /// @return the KeyVal object corresponding to \c path ;
  ///         if \c path does not exist, return an empty KeyVal
  KeyVal keyval(const key_type& path) const {
    auto abs_path = resolve_path(path);
    return KeyVal(top_tree_, class_registry_, abs_path);
  }

  /// \brief returns a shared_ptr to the (top) tree
  std::shared_ptr<ptree> top_tree() const { return top_tree_; }
  /// \brief returns a shared_ptr to this (sub)tree
  /// @param path the path to the subtree
  /// @note the result is aliased against the top tree
  std::shared_ptr<ptree> tree() const;

  /// checks whether the given path exists
  /// @param path the path
  /// @return true if \c path exists
  bool exists(const key_type& path) const {
    return exists_(resolve_path(path));
  }

  /// check whether the given class exists
  /// @param path the path
  /// @return true if \c path class exists
  bool exists_class(const key_type& path) const{
    bool exist_class = false;
    auto cptr = class_registry_->find(resolve_path(path));
    if (cptr != class_registry_->end()){
      exist_class = true;
    }
    return exist_class;
  }

  /// counts the number of children of the node at this path
  /// @param path the path
  /// @return 0 if \c path does not exist or it points to a simple keyword,
  ///         the number of elements (>=0) if \c path points to an array or
  ///         a keyword group
  size_t count(const key_type& path) const {
    auto resolved_path = resolve_path(path);
    auto child_opt = top_tree_->get_child_optional(ptree::path_type{path, separator});
    if (child_opt == boost::optional<ptree&>())
      return 0;
    else
      return child_opt->size();
  }

  /// assign simple \c value at the given path (overwrite, if necessary)
  /// @param value a "simple" value, i.e. it can be converted to
  /// a KeyVal::key_type using a
  /// std::basic_ostream<KeyVal::key_type::value_type>
  template <typename T,
            typename = enable_if_t<not KeyVal::is_sequence<T>::value >>
  KeyVal& assign(const key_type& path, const T& value) {
    auto abs_path = to_absolute_path(path);
    top_tree_->put(ptree::path_type{abs_path, separator}, value);
    return *this;
  }

  /// assign a sequence container at the given path (overwrite, if necessary)
  /// @tparam SequenceContainer any container for which
  /// KeyVal::is_sequence<SequenceContainer> is a \c std::true_type,
  ///         currently any of the following is allowed: \c std::array, \c
  ///         std::vector, \c std::list .
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
      enable_if_t<KeyVal::is_sequence<SequenceContainer>::value>* = nullptr) {
    auto abs_path = to_absolute_path(path);
    ptree obj;
    size_t count = 0;
    for (const auto& v : value) {
      ptree pt;
      auto key =
          json_style
              ? key_type("")
              : std::to_string(count);  // assumes key_type == std::string
      pt.put("", v);
      obj.push_back(std::make_pair(key, pt));
      ++count;
    }
    top_tree_->add_child(ptree::path_type{abs_path, separator}, obj);
    return *this;
  }

  /// assign the given pointer to a DescribedClass at the given path (overwrite,
  /// if necessary)
  /// @warning these key/value pairs are not part of ptree, hence cannot be
  /// written to JSON/XML
  KeyVal& assign(const key_type& path,
                 const std::shared_ptr<DescribedClass>& value) {
    auto abs_path = to_absolute_path(path);
    (*class_registry_)[abs_path] = value;
    return *this;
  }

  /// return the result of converting the value at the given path to the desired
  /// type.

  /// @tparam T the desired value type, must be a "simple" type that can be
  /// accepted by KeyVal::assign()
  /// @param path the path to the value
  /// @return value stored at \c path converted to type \c T
  /// @throws KeyVal::bad_input if path not found or cannot convert value
  /// representation to the desired type
  template <typename T,
            typename = enable_if_t<not KeyVal::is_sequence<T>::value>>
  T value(const key_type& path) const {
    auto abs_path = resolve_path(path);
    T result;
    try {
      result = top_tree_->get<T>(ptree::path_type{abs_path, separator});
    } catch (boost::property_tree::ptree_bad_data& x) {
      throw KeyVal::bad_input(
          std::string("value ") + x.data<KeyVal::value_type>() +
              " could not be converted to the desired datatype",
          abs_path);
    }
    return result;
  }

  /// return value corresponding to a path and convert to a std::vector.
  /// @tparam T the desired value type
  /// @param path the path
  /// @throw KeyVal::bad_input if \c path is bad
  /// @return value of type \c T
  template <typename SequenceContainer>
  SequenceContainer value(
      const key_type& path,
      enable_if_t<KeyVal::is_sequence<SequenceContainer>::value>* =
          nullptr) const {
    auto abs_path = resolve_path(path);
    using value_type = typename SequenceContainer::value_type;
    SequenceContainer result;
    try {
      auto vec_ptree =
          top_tree_->get_child(ptree::path_type{abs_path, separator});
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
        *iter = elem_ptree.second.get_value<value_type>();
        ++iter;
        ++count;
      }
      return result;
    } catch (boost::property_tree::ptree_bad_data& x) {
      throw KeyVal::bad_input(
          std::string("value ") + x.data<KeyVal::value_type>() +
              " could not be converted to the desired datatype",
          abs_path);
    } catch (boost::property_tree::ptree_bad_path& x) {
      throw KeyVal::bad_input(
          std::string("path ") + x.path<key_type>() + " not found", abs_path);
    }
  }

  /// return value corresponding to a path and convert to the desired type.
  /// @tparam T the desired value type
  /// @param path the path
  /// @param default_value
  /// @return value of type \c T
  template <typename T,
            typename = enable_if_t<not KeyVal::is_sequence<T>::value>>
  T value(const key_type& path, const T& default_value) const {
    auto abs_path = resolve_path(path);
    T result;
    if (not exists_(abs_path))
      result = default_value;
    else
      result = top_tree_->get<T>(ptree::path_type{abs_path, separator});
    return result;
  }

  /// return a pointer to an object at the given path

  /// Constructs an object of the type given by the value of
  /// keyword "type" in this KeyVal object. In order to be constructible
  /// using this method a type must be derived from DescribedClass <i>and</i>
  /// registered.
  /// @tparam T a class derived from DescribedClass corresponding to the value
  /// of path:"type" or its base
  /// @param path the path
  /// @return std::shared_Ptr \c T
  /// @throws KeyVal::bad_input if key not found or cannot convert value
  /// representation to the desired type
  template <typename T, typename = enable_if_t<Describable<T>::value>>
  std::shared_ptr<T> class_ptr(const key_type& path = key_type()) const {
    // if this class already exists in the registry under path
    // (e.g. if the ptr was assigned programmatically), return immediately
    {
      auto cptr = class_registry_->find(path);
      if (cptr != class_registry_->end())
        return std::dynamic_pointer_cast<T>(cptr->second);
    }

    // otherwise, resolve the path and build the class (or pick up the cached
    // copy)
    auto abs_path = resolve_path(path);

    // if this class already exists, return the ptr
    auto cptr = class_registry_->find(abs_path);
    if (cptr != class_registry_->end())
      return std::dynamic_pointer_cast<T>(cptr->second);

    // get class name
    std::string derived_class_name;
    auto type_path = concat_path(abs_path, "type");
    if (not exists_(type_path))
      throw KeyVal::bad_input(
          std::string("missing \"type\" in object at path ") + path, abs_path);
    derived_class_name =
        top_tree_->get<std::string>(ptree::path_type{type_path, separator});

    // map class type to its ctor
    auto ctor = DescribedClass::type_to_keyval_ctor(derived_class_name);
    auto result = (*ctor)(*this);
    (*class_registry_)[abs_path] = result;
    return std::dynamic_pointer_cast<T>(result);
  }

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

 private:
  std::shared_ptr<ptree> top_tree_;
  using class_registry_type =
      std::map<std::string, std::shared_ptr<DescribedClass>>;
  std::shared_ptr<class_registry_type> class_registry_;
  key_type path_;  //!< path from the top of \c top_tree_ to this subtree

  /// creates a KeyVal
  KeyVal(const std::shared_ptr<ptree>& top_tree,
         const std::shared_ptr<class_registry_type>& class_registry,
         const key_type& path)
      : top_tree_(top_tree), class_registry_(class_registry), path_(path) {}

  /// creates a KeyVal using am existing \c ptree
  /// @param pt a ptree
  KeyVal(ptree pt)
      : top_tree_(std::make_shared<ptree>(std::move(pt))),
        class_registry_(),
        path_("") {}

  /// given a path that contains ".." elements, returns the equivalent path
  /// without such elements
  /// @example realpath("tmp:..:x") returns "x"
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

  /// @returns true if \c path did not include a reference, even if path does not exist
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

  /// checks whether the given path exists; does not resolve path unlike
  /// KeyVal::exists()
  /// @param path the path
  /// @return true if \c path exists
  bool exists_(const key_type& path) const {
    return top_tree_->get_child_optional(ptree::path_type{path, separator}) !=
           boost::optional<ptree&>();
  }

 public:
  /// \brief KeyVal exception
  struct bad_input : public std::runtime_error {
    bad_input(const std::string& _what, const key_type& path)
        : std::runtime_error(_what + "(path=" + path + ")") {}
    virtual ~bad_input() noexcept {}
  };

};  // KeyVal

/// union of two KeyVal objects
/// @note obsoletes sc::AggregateKeyVal
/// @ingroup KeyValCore
KeyVal operator+(const KeyVal& first, const KeyVal& second);
}

/// @}

#endif /* SRC_MPQC_UTIL_KEYVAL_KEYVAL_HPP_ */
