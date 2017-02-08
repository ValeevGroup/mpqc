#ifndef SRC_MPQC_UTIL_MISC_PROVIDER_H_
#define SRC_MPQC_UTIL_MISC_PROVIDER_H_

#include <memory>

// This provides syntactic sugar for Acyclic Visitor.
// We are not using the traditional names since they only confuse users.
// Translation: Provider = Visitor, Property = Visitable

namespace mpqc {

/// \brief Base for classes that provide \c Properties .

/// This provides to the class that inherits this an ability to visit
/// each property \c P in \c Properties by overloading
/// the corresponding \c P::Provider::can_evaluate and \c
/// P::Provider::evaluate methods.
/// @tparam Properties the property type list
template <typename... Properties>
class Provides : public Properties::Provider... {};

/// @return true if Provider can provide Property
template <typename Property, typename Provider>
bool provides(const std::shared_ptr<Provider>& provider) {
  return std::dynamic_pointer_cast<Property::Provider>(provider) == nullptr;
}

namespace detail {
/// has_provider<T>::value is true if T::Provider is a valid type
template <typename T>
struct has_provider {
  template <typename U>
  static std::true_type test(typename U::Provider*);
  template <typename U>
  static std::false_type test(...);
  static constexpr const bool value = decltype(test<T>(nullptr))::value;
};

/// to_pointer<T> converts obj to a pointer:
/// - if T is a pointer type, will return obj
/// - if T is a value or a reference, will take the address
/// - if T is a shared_ptr, will return the corresponding raw pointer
template <typename T>
T* to_pointer(T* obj) {
  return obj;
}
template <typename T>
auto to_pointer(T& obj) -> std::enable_if_t<
    utility::meta::is_shared_ptr<typename std::decay<T>::type>::value,
    decltype(obj.get())> {
  return obj.get();
}
template <typename T>
std::enable_if_t<
    !utility::meta::is_shared_ptr<typename std::decay<T>::type>::value,
    typename std::decay<T>::type*>
to_pointer(T& obj) {
  return &obj;
}

/// obtains a description of the object pointer to by \c obj_ptr
template <typename T>
std::string description(T* obj_ptr) {
  auto dc_obj_ptr = dynamic_cast<DescribedClass*>(obj_ptr);
  if (dc_obj_ptr != nullptr) {
    return std::string("class ") + dc_obj_ptr->class_key();
  } else {
    std::ostringstream oss;
    oss << "object @ " << obj_ptr << " of type " << typeid(*obj_ptr).name();
    return oss.str();
  }
}

}  // namespace detail

/// Evaluates \c property using \c provider
template <typename Property, typename Provider, typename... EvaluateArgs>
std::enable_if_t<(detail::has_provider<Property>::value &&
                  !std::is_const<Property>::value),
                 void>
evaluate(Property& property, Provider& provider, EvaluateArgs... eval_args) {
  auto provider_ptr = detail::to_pointer(provider);
  auto* evaluator = dynamic_cast<typename Property::Provider*>(provider_ptr);
  if (evaluator == nullptr) {
    std::ostringstream oss;
    oss << detail::description(provider_ptr) << " does not compute "
        << detail::description(&property)
        << ", needs to derive from an appropriate Provides<>";
    throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
  evaluator->evaluate(&property, eval_args...);
}

/// Evaluates \c property using \c provider
template <typename Property, typename Provider>
std::enable_if_t<detail::has_provider<Property>::value, Property&> operator<<(
    Property& property, Provider& provider) {
  evaluate(property, provider);
  return property;
}

}  // namespace mpqc

#endif /* SRC_MPQC_UTIL_MISC_PROVIDER_H_ */
