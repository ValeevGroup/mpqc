#ifndef SRC_MPQC_UTIL_MISC_OBSERVER_H_
#define SRC_MPQC_UTIL_MISC_OBSERVER_H_

#include <map>
#include <vector>

namespace mpqc {
namespace utility {

/// @brief helps to set up messaging between objects via the Observer pattern

/// Use as a base class for any class that wants to send messages to other
/// recipients.
/// Recipients register callbacks by calling Observable::register_message .
/// @note Observable holds weak (non-owning) pointers to the observers. Message
/// lifetime is managed automatically as long as the observers derive from the
/// Observer class. Observer will hold std::shared_ptr to this, hence the need for
/// the std::enable_shared_from_this base.
/// @sa Observer
/// @tparam Derived the derived class for which std::shared_ptr exists (as required
///                 for std::enable_shared_from_this to be usable)
template <typename Derived>
class Observable : public std::enable_shared_from_this<Derived> {
 public:
  typedef std::map<void *, std::function<void()>> messages_type;

  /// adds a message to the list, returns a receipt that the observer
  /// needs to destroy at the end of its life (i.e., in the dtor).
  template <typename Observer>
  std::shared_ptr<void *const> register_message(Observer *observer,
                                                std::function<void()> message) {
    void *observer_void_ptr = static_cast<void *>(observer);
    messages_.insert(std::make_pair(static_cast<void *>(observer), message));
    auto ptr_to_this = this->shared_from_this();  // TODO move-capture ptr_to_this in C++14
    std::shared_ptr<void *> result(new void *, [ptr_to_this](void **ptr) {
      ptr_to_this->messages_.erase(*ptr);
      delete ptr;
    });
    *result = observer_void_ptr;
    return result;
  }

 protected:
  void message() {
    for (auto &message : messages_) {
      (message.second)();
    }
  }

 private:
  // bag of messages to call every time Observable::message() is called
  messages_type messages_;
};

/// @brief helps to set up messaging between objects via the Observer pattern

/// Use as a base class for any class that wants to receive messages from
/// Observable objects.
/// @sa Observable
class Observer {
 public:
  /// The default ctor
  /// @note using the default ctor requires you to call register_message
  /// manually
  Observer() = default;

  /// This constructor registers one message
  /// @param observee the Observee object to send us the \c message
  /// @param message the message to be sent
  template <typename Observee>
  Observer(Observee *observee, std::function<void()> message) {
    register_message(observee, message);
  }

  // cannot copy, must construct a "copy" manually
  Observer(const Observer &other) = delete;
  Observer &operator=(const Observer &other) = delete;

  Observer(Observer &&other) = default;
  Observer &operator=(Observer &&other) = default;

  virtual ~Observer() { }

 protected:
  /// registers a message with an Observee object, will be called whenever
  /// Observee::update()
  /// is called on \c observee
  /// @tparam Observee a class derived from Observable
  template <typename Observee>
  void register_message(Observee *observee, std::function<void()> message) {
    receipts_.emplace_back(std::move(observee->register_message(this, message)));
  }

  void clear_messages() { receipts_.clear(); }

 private:
  std::vector<std::shared_ptr<void *const>> receipts_;
};

}  // namespace utility
}  // namespace mpqc

#endif /* SRC_MPQC_UTIL_MISC_OBSERVER_H_ */
