/*
 * keyval.hpp
 *
 *  Created on: Apr 18, 2016
 *      Author: evaleev
 */

#ifndef SRC_MPQC_UTIL_KEYVAL_KEYVAL_HPP_
#define SRC_MPQC_UTIL_KEYVAL_KEYVAL_HPP_

#include <map>
#include <string>
#include <memory>
#include <type_traits>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace mpqc {

  namespace exception {
    struct bad_input : public std::runtime_error {
        bad_input(const char* _what) : std::runtime_error(_what) {}
        bad_input(const std::string& _what) : std::runtime_error(_what) {}
    };
  }

  class KeyVal;

  /// This class helps construction of (smart pointers to) C++ objects from , e.g., text input

  /// To construct a class you want to be able to construct from KeyVal objects do this:
  /// 1. make DescribedClass a public base for this class, and 2. create a registrar object
  /// for the class as a class static member or any other object type with static storge duration,
  /// e.g. a variable declared in global namespace.
  class DescribedClass {
    public:
      DescribedClass() = default;
      virtual ~DescribedClass() {}

      typedef std::shared_ptr<DescribedClass> (*keyval_ctor_wrapper_type)(const KeyVal&);
      static keyval_ctor_wrapper_type type_to_keyval_ctor(const std::string& type_name) {
        auto& registry = keyval_ctor_registry();
        if (registry.find(type_name) == registry.end())
          throw std::runtime_error("DescribedClass::type_to_keyval_ctor -- type not registered");
        return registry[type_name];
      }
      template <typename T> static void register_keyval_ctor(const std::string& type_name) {
        auto& registry = keyval_ctor_registry();
        assert(registry.find(type_name) == registry.end());
        registry[type_name] = keyval_ctor_wrapper<T>;
      }

      /// make a single instance of this to register the KeyVal ctor of type \c T
      template <typename T> struct registrar {
          registrar(const std::string& type_name) {
            DescribedClass::register_keyval_ctor<T>(type_name);
          }
      };

    private:
      using keyval_ctor_registry_type = std::map<std::string,keyval_ctor_wrapper_type>;
      // this is needed to force registry initialization BEFORE its use
      static keyval_ctor_registry_type& keyval_ctor_registry() {
        static keyval_ctor_registry_type keyval_ctor_registry_;
        return keyval_ctor_registry_;
      }
      template <typename T> static std::shared_ptr<DescribedClass> keyval_ctor_wrapper(const KeyVal& kv) {
        return std::make_shared<T>(kv);
      }
  };

  /**
   * @ingroup CoreKeyVal
   The KeyVal class is designed to simplify the process of allowing
   a user to specify keyword/value associations to a C++ program.  A
   flexible input style and ease of use for the programmer is achieved with
   this method.  Keywords are represented by null terminated character arrays.
   The keywords are organized hierarchically, in a manner similar to the way
   that many file systems are organized.  One character is special,
   ":", which is used to separate the various hierarchical labels,
   which are referred to as "segments", in the keyword.

   A convention for specifying arrays is provided by KeyVal.  Each
   index of the array is given by appending a segment containing the
   character representation of the index.  Thus, "array:3:4" would be
   a the keyword corresponding to fourth row and fifth column of
   "array", since indexing starts at zero.

   To allow the KeyVal class to have associations that can represent
   data for classes, the keyword can be associated with a class as well as
   a value.  This permits polymorphic data to be unambiguously represented
   by keyword/value associations.  Most use of KeyVal need not be
   concerned with this.
  */

  ///  a (sub)tree of Key=Value pairs implemented as a Boost PropertyTree.

  /// KeyVal is (a subtree of) a tree of key-value pairs.
  ///
  /// @note Since this is default-constructible and directly mutable,
  ///       this obsoletes sc::AssignedKeyVal
  class KeyVal {
    public:
      /// data type for representing a property tree
      using ptree = boost::property_tree::iptree;
      using key_type = ptree::key_type; // = std::string
      constexpr static char separator = ':';

      /// creates empty KeyVal
      KeyVal() : top_tree_(std::make_shared<ptree>()),
          class_registry_(std::make_shared<class_registry_type>()),
          path_("") {}

      /// constructs from a JSON file
      /// @note obsoletes sc::ParsedKeyVal
      KeyVal(std::string filename);

      /// (mostly shallow) copy ctor (to clone top-level tree use top_tree())
      KeyVal(const KeyVal& other) = default;
      /// copy assignment
      KeyVal& operator=(const KeyVal& other) = default;

      /// extract a subtree KeyVal located at \c path
      /// @note obsoletes PrefixKeyVal
      /// @note the sub-tree should not contain references
      /// @param path the path to the subtree; absolute (within the top_tree) and
      ///        relative (with respect to this subtree) paths are allowed.
      /// @return the KeyVal object corresponding to \c path ;
      ///         if \c path does not exist, return an empty KeyVal
      KeyVal keyval(const key_type& path) const {
        auto abs_path = resolve_path(path);
        return KeyVal(top_tree_, class_registry_, abs_path);
      }

      /// returns a shared_ptr to the top tree
      std::shared_ptr<ptree> top_tree() const {
        return top_tree_;
      }
      /// returns a shared_ptr to this (sub)tree
      /// @param path the path to the subtree
      /// @note the result is aliased against the top tree
      std::shared_ptr<ptree> tree() const {
        std::shared_ptr<ptree> result(top_tree_, &top_tree_->get_child(ptree::path_type{path_, separator}));
        return result;
      }

      /// assign simple \c value to key with path \c path ; \c value is simple if it can be converted to
      /// a \c std::string using a std::ostream
      template <typename T>
      KeyVal& assign(const key_type& path, const T& value) {
        auto abs_path = to_absolute_path(path);
        top_tree_->put(ptree::path_type{abs_path, separator}, value);
        return *this;
      }

      /// assign a \c std::vector to key with path \c path
      template <typename T>
      KeyVal& assign(const std::string& path, const std::vector<T>& value) {
        auto abs_path = to_absolute_path(path);
        ptree obj;
        for(const auto& v: value) {
          ptree pt;
          pt.put("", v);
          obj.push_back(std::make_pair("", pt));
        }
        top_tree_->add_child(ptree::path_type{abs_path, separator}, obj);
        return *this;
      }

      /// assign a ptr to a DescribedClass to key with path \c path
      KeyVal& assign(const std::string& path, const std::shared_ptr<DescribedClass>& value) {
        assert(false); // not yet implemented
      }

      /// return value corresponding to a path and convert to the desired type.
      /// @tparam T the desired value type
      /// @param path the key
      /// @return value of type \c T
      /// @throws mpqc::exception::bad_input if path not found or cannot convert value representation to the desired type
      template <typename T>
      T value(const key_type& path) const {
        auto abs_path = resolve_path(path);
        T result;
        try {
          result = top_tree_->get<T>(ptree::path_type{abs_path, separator});
        }
        catch(boost::property_tree::ptree_bad_data&) {
          throw mpqc::exception::bad_input(abs_path);
        }
        return result;
      }

      /// return value corresponding to a path and convert to the desired type.
      /// @tparam T the desired value type
      /// @param path the path
      /// @param default_value
      /// @return value of type \c T
      template <typename T>
      T value(const key_type& path, const T& default_value) const {
        auto abs_path = resolve_path(path);
        T result;
        try {
          result = top_tree_->get<T>(ptree::path_type{abs_path, separator});
        }
        catch(boost::property_tree::ptree_bad_data&) {
          result = default_value;
        }
        return result;
      }

      /// return a pointer to an object specified by this KeyVal

      /// Constructs an object of the type given by the value of
      /// keyword "type" in this KeyVal object. In order to be constructible
      /// using this method a type must be derived from DescribedClass and
      /// it must be registered by calling
      /// DescribedClass::register_keyval_ctor() with the type name
      /// given as the string argument.
      /// @tparam T the class corresponding to the "type" keyword
      ///           or its base
      /// @param path the path
      /// @return value of type \c T
      /// @throws mpqc::exception::bad_input if key not found or cannot convert value representation to the desired type
      template <typename T,
                typename = typename std::enable_if<std::is_base_of<DescribedClass,T>::value>::type>
      std::shared_ptr<T> class_ptr() const {

        key_type abs_path = path_;

        // if this class already exists, return the ptr
        auto cptr = class_registry_->find(abs_path);
        if (cptr != class_registry_->end())
          return std::dynamic_pointer_cast<T>(cptr->second);

        // get class name
        std::string derived_class_name;
        try {
          auto type_path = concat_path(abs_path,"type");
          derived_class_name = top_tree_->get<std::string>(ptree::path_type{type_path, separator});
        }
        catch(boost::property_tree::ptree_bad_data&) {
          throw mpqc::exception::bad_input(std::string("KeyVal: missing \"type\" in object at "));
        }

        // map class type to its ctor
        auto ctor = DescribedClass::type_to_keyval_ctor(derived_class_name);
        auto result = (*ctor)(*this);
        (*class_registry_)[abs_path] = result;
        return std::dynamic_pointer_cast<T>(result);
      }

      /// write to stream in XML form
      void
      write_xml(std::basic_ostream< typename key_type::value_type > & os) const {
        boost::property_tree::xml_parser::write_xml(os, *(tree()));
      }
      /// write from stream in XML form
      void
      read_xml(std::basic_istream< typename key_type::value_type > & is) const {
        boost::property_tree::xml_parser::read_xml(is, *(tree()));
      }

      /// write to stream in JSON form
      void
      write_json(std::basic_ostream< typename key_type::value_type > & os) const {
        boost::property_tree::json_parser::write_json(os, *(tree()));
      }
      /// read from stream in JSON form
      void
      read_json(std::basic_istream< typename key_type::value_type > & is) const {
        boost::property_tree::json_parser::read_json(is, *(tree()));
      }

    private:
      std::shared_ptr<ptree> top_tree_;
      using class_registry_type = std::map<std::string,std::shared_ptr<DescribedClass>>;
      std::shared_ptr<class_registry_type> class_registry_;
      key_type path_;  //!< path from the top of \c top_tree_ to this subtree

      /// creates a KeyVal
      KeyVal(const std::shared_ptr<ptree>& top_tree,
             const std::shared_ptr<class_registry_type>& class_registry,
             const key_type& path) :
        top_tree_(top_tree),
        class_registry_(class_registry),
        path_(path) {}

      /// creates a KeyVal using am existing \c ptree
      /// @param pt a ptree
      KeyVal(ptree pt) :
        top_tree_(std::make_shared<ptree>(std::move(pt))),
        class_registry_(),
        path_("") {}

//      /// This replaces values that are references with the actual values
//      /// thus converting a DAG into a tree
//      /// @note repeatedly recurses, hence can be expensive
//      /// @note circular references are not detected, thus to prevent infinite recursion possible
//      ///       will throw if too many iterations are needed
//      /// @param niter_max max number of iterations
//      /// @throws \c mpqc::exception::bad_input if max number of iterations exceeded
//      void resolve_refs(size_t niter_max = 10) {
//        size_t niter = 0;
//        while(niter<niter_max) {
//          bool done = not has_refs();
//          if (done)
//            return;
//          resolve_refs_once();
//          ++niter;
//        }
//        // too many iterations if still here
//        throw mpqc::exception::bad_input(std::string("KeyVal: excessive or circular references"));
//      }
//
//      bool has_refs_helper(size_t depth, const ptree& pt) {
//        if (pt.size() > 0) { // non-leaf
//          for(auto& kv_pair: pt) {
//            ptree subtree = kv_pair.second;
//            auto result = has_refs_helper(depth+1, subtree);
//            if (result) return true;
//          }
//        }
//        else // leaf
//          return pt.data().size() != 0 && pt.data()[0] == '$';
//        return false;
//      }
//      /// @return true if any value is a ref, i.e. starts with a "$"
//      bool has_refs() {
//        return has_refs_helper(0, pt_);
//      }
//
//      void resolve_refs_once_helper(size_t depth, ptree& top_pt, ptree& pt) {
//        if (pt.size() > 0) { // non-leaf
//          for(auto& kv_pair: pt) {
//            ptree& subtree = kv_pair.second;
//            resolve_refs_once_helper(depth+1, top_pt, subtree);
//          }
//        }
//        else { // leaf
//          auto value_str = pt.data();
//          auto is_a_ref = value_str.size() != 0 && value_str[0] == '$';
//          if (is_a_ref) {
//            // only can handle absolute references that start with "$:"
//            assert(value_str[1] == separator);
//            // assuming no ".."
//            auto abs_ref_key = value_str.substr(2); //drop the leading "$:"
//            ptree& ref_pt = top_pt.get_child(ptree::path_type{abs_ref_key, separator});
//            pt = ref_pt;
//          }
//        }
//      }
//      void resolve_refs_once() {
//        resolve_refs_once_helper(0, pt_, pt_);
//      }

      /// given a path that contains ".." elements, returns the equivalent path without such elements
      /// @example realpath("tmp:..:x") returns "x"
      /// @throw mpqc::exception::bad_input if path is invalid; an example is ".." .
      static key_type realpath(const key_type& path) {
        key_type result = path;
        const key_type parent_token("..");
        auto parent_token_location = result.find(parent_token);
        while (parent_token_location != key_type::npos) {
          // make sure .. is preceeded by :
          if (parent_token_location == 0 || path[parent_token_location-1] != separator) {
            throw mpqc::exception::bad_input("KeyVal::realpath() -- bad path");
          }
          // did we find '..:' or '..\0'?
          auto parent_token_length = (result.find(separator,parent_token_location+1) == parent_token_location+2) ? 3 : 2;
          // find the beginning of the preceding path element
          // NB handle the corner case of path starting with ":.."
          auto preceding_token_location = parent_token_location > 1 ? path.rfind(separator, parent_token_location-2) : 0;
          if (preceding_token_location == key_type::npos) // if separator not found, the preceeding element is the beginning of the path
            preceding_token_location = 0;
          result.erase(preceding_token_location, (parent_token_location + parent_token_length - preceding_token_location));
          parent_token_location = result.find(parent_token);
        }
        return result;
      }

      /// computes an absolute path from a given path prefix and an (absolute or relative) path
      /// @param path_prefix
      /// @param path an absolute path (starts with "$:") or a relative path (the leading "$" is dropped)
      /// @return the absolute path
      /// @note this does not resolve references
      static key_type to_absolute_path(const key_type& path_prefix, const key_type& path) {
        auto is_ref = path.size() != 0 && path[0] == '$';
        auto is_abs_ref = is_ref && path.size() > 1 && path[1] == ':';
        key_type abs_path;
        key_type path_copy = path;
        if (is_abs_ref) {
          abs_path = path_copy.substr(2); //drop the leading "$:"
        }
        else {
          key_type rel_path = is_ref ? path_copy.substr(1) : path_copy; //drop the leading "$"
          abs_path = concat_path(path_prefix,rel_path);
        }
        abs_path = realpath(abs_path);
        return abs_path;
      }

      /// computes an absolute path from a given (absolute or relative) path
      /// @param path an absolute path (starts with "$:") or a relative path (the leading "$" is dropped)
      /// @return the absolute path
      /// @note this does not resolve references
      key_type to_absolute_path(const key_type& path) const {
        return to_absolute_path(path_, path);
      }

      /// Concatenates two paths
      /// @return a path obtained by concatenating \c prefix and \c postfix ;
      static key_type concat_path(const key_type& prefix, const key_type& postfix) {
        //if (postfix.size() > 0) assert(postfix[0] != separator);
        auto result = (prefix == "") ? postfix : prefix + separator + postfix;
        return result;
      }

      /// normalizes path by 1. converting path to absolute path (see \c to_absolute_path())
      /// and 2. resolving any refs in the path
      key_type resolve_path(const key_type& path) const {
        auto abs_path = to_absolute_path(path);
        auto result = resolve_refs(abs_path);
        return result;
      }

      /// This resolves the references in the path.
      /// @note repeatedly recurses, hence can be expensive
      /// @note circular references are not detected, thus to prevent infinite recursion possible
      ///       will throw if too many iterations are needed
      /// @param niter_max max number of iterations
      /// @throws \c mpqc::exception::bad_input if max number of iterations exceeded
      key_type
      resolve_refs (const key_type& path, size_t niter_max = 10) const {
        key_type result = path;
        size_t niter = 0;
        while (niter < niter_max) {
          bool done = resolve_first_ref(result);
          if (done) return result;
          ++niter;
        }
        // too many iterations if still here
        throw mpqc::exception::bad_input (
            std::string ("KeyVal: excessive or circular references"));
      }

      /// @returns true if \c path did not include a reference
      bool resolve_first_ref(key_type& path) const {
        if (path.size() == 0) return true; // handle this corner case now to make the rest a bit cleaner

        auto subpath_last_char = path.find(separator, 1); // if path starts ":" no need to check root of the tree
        if (subpath_last_char == key_type::npos) subpath_last_char = path.size(); // this ensures that the full path is checked
        while(subpath_last_char <= path.size()) {
          auto subpath = path.substr(0,subpath_last_char);
          auto value = top_tree_->get<key_type>(ptree::path_type{subpath, separator}); // value corresponding to subpath
          auto subpath_is_a_ref = (value.size() > 0 && value[0] == '$');
          if (subpath_is_a_ref) {
            // resolve subpath key + value to a new subpath
            auto subpath_base = path_pop(subpath);
            auto new_subpath = to_absolute_path(subpath_base, value);

            // replace the old subpath with the new subpath
            path = new_subpath + path.substr(subpath.size());

            return false;
          }
          if (subpath_last_char == path.size()) // subpath = path and it's not a ref? done
            return true;
          subpath_last_char = path.find(separator, subpath_last_char+1);
          if (subpath_last_char == key_type::npos) subpath_last_char = path.size(); // this ensures that the full path is checked
        }
        return true; // no refs found
      }

      /// pops off last ":token" (if separator not found, the result is an empty path)
      static key_type path_pop(const key_type& path) {
        key_type result = path;
        auto last_separator_location = result.rfind(separator);
        if (last_separator_location == key_type::npos) last_separator_location = 0;
        result.erase(last_separator_location);
        return result;
      }

  }; // KeyVal

  /// union of two KeyVal objects
  /// @note obsoletes sc::AggregateKeyVal
  KeyVal operator+(const KeyVal& first, const KeyVal& second);


}

#endif /* SRC_MPQC_UTIL_KEYVAL_KEYVAL_HPP_ */
