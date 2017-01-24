#ifndef SRC_MPQC_UTIL_MISC_OBSERVER_H_
#define SRC_MPQC_UTIL_MISC_OBSERVER_H_

#include <map>
#include <vector>

namespace mpqc {
namespace utility {

/// @brief helps to set up messaging between objects via the Observer pattern

/// Use as a base class for any class that wants to send messages to other recipients.
/// Recipients register callbacks by calling Observable::register_message .
/// @note Observable holds weak (non-owning) pointers to the observers. Message lifetime
///       is managed automatically as long as the observers derive from the Observer class.
/// @sa Observer
class Observable {
 public:
  typedef std::map<void *, std::function<void()>> messages_type;

  /// adds a message to the list, returns a receipt that the observer
  /// needs to destroy at the end of its life (i.e., in the dtor).
  template <typename Observer>
  std::shared_ptr<void * const> register_message(Observer *observer,
                                                 std::function<void()> message) {
    void *observer_void_ptr = static_cast<void *>(observer);
    messages_.insert(std::make_pair(static_cast<void *>(observer), message));
    std::shared_ptr<void *> result(new void *, [this](void **ptr) {
      messages_.erase(*ptr);
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

/// Use as a base class for any class that wants to receive messages from Observable objects.
/// @sa Observable
class Observer {
 public:
  /// The default ctor
  /// @note using the default ctor requires you to call register_message manually
  Observer() = default;

  /// This constructor registers one message
  /// @tparam Observee a class derived from Observable
  /// @param observee the Observee object to send us the \c message
  /// @param message the message to be sent
  template <typename Observee>
  Observer(Observee *observee, std::function<void()> message) {
    register_message(observee, message);
  }

  virtual ~Observer() = default;

 protected:
  /// registers a message with a Messenger object, will be called whenever MessengerBase::update()
  /// is called on \c messenger
  /// @tparam Observee a class derived from Observable
  template <typename Observee>
  void register_message(Observee *observee, std::function<void()> message) {
    receipts_.emplace_back(std::move(observee->register_message(this, message)));
  }

  void clear_messages() { receipts_.clear(); }

 private:
  std::vector<std::shared_ptr<void * const>> receipts_;
};

}  // namespace utility
}  // namespace mpqc



#endif /* SRC_MPQC_UTIL_MISC_OBSERVER_H_ */
