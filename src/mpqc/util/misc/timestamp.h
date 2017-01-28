#ifndef SRC_MPQC_UTIL_MISC_TIMESTAMP_H_
#define SRC_MPQC_UTIL_MISC_TIMESTAMP_H_

#include <atomic>
#include <iostream>
#include <memory>

namespace mpqc {
namespace utility {

/// \brief Produces timestamps

/// Need to recompute is detected by comparing timestamps of property parameters
/// and value
/// timestamps are produced by this factory method
class TimestampFactory {
 public:
  using timestamp_type = uint64_t;
  static timestamp_type make() {
    static std::atomic<timestamp_type> current_timestamp{0};
    return current_timestamp++;
  }
};

/// Timestampable<T> is a proxy to T that keeps the timestamp of the last
/// modification.
/// Stores T on heap to support the default_initialized state.
template <typename T>
class Timestampable {
 public:
  using timestamp_type = TimestampFactory::timestamp_type;

  Timestampable() : timestamp_(get_timestamp()) {}
  explicit Timestampable(T&& val)
      : value_(std::make_shared<T>(val)), timestamp_(get_timestamp()) {}
  explicit Timestampable(std::shared_ptr<T> val)
      : value_(val), timestamp_(get_timestamp()) {}

  /// copy ctor keeps the timestamp, deep copies value
  Timestampable(const Timestampable& other)
      : timestamp_(other.timestamp_), value_(std::make_shared<T>(T(other))) {}

  Timestampable& operator=(const T& val) {
    value_ = std::make_shared<T>(val);
    timestamp_ = get_timestamp();
    return *this;
  }
  Timestampable& operator=(std::shared_ptr<T> val) {
    value_ = val;
    timestamp_ = get_timestamp();
    return *this;
  }

  /// retrieve a (non-const) reference value_ updates the timestamp
  operator T&() {
    assert(value_ != nullptr);
    timestamp_ = get_timestamp();
    return *value_;
  }

  /// retrieve a const reference value_
  operator const T&() const {
    assert(value_ != nullptr);
    return *value_;
  }

  /// retrieve a shared pointer, updates the timestamp
  operator std::shared_ptr<T>() {
    assert(value_ != nullptr);
    timestamp_ = get_timestamp();
    return value_;
  }

  /// retrieve a const shared pointer
  operator std::shared_ptr<const T>() const {
    assert(value_ != nullptr);
    return value_;
  }

  /// @return the current timestamp
  const timestamp_type& timestamp() const { return timestamp_; }
  /// @return the current value
  std::shared_ptr<const T> value() const { return value_; }

 private:
  std::shared_ptr<T> value_;
  timestamp_type timestamp_;

  static timestamp_type get_timestamp() { return TimestampFactory::make(); }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Timestampable<T>& x) {
  os << static_cast<const T&>(x);
  return os;
}

}  // namespace utility
}  // namespace mpqc

#endif /* SRC_MPQC_UTIL_MISC_TIMESTAMP_H_ */
