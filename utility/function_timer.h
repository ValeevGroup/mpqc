#pragma once
#ifndef TCC_UTILITY_FUNCTIONTIMER_H
#define TCC_UTILITY_FUNCTIONTIMER_H

#include <chrono>
#include <utility>

namespace tcc {
namespace utility {

template <typename Fn>
class FunctionTimer {
  public:
    FunctionTimer(Fn fn) : fn_{fn} {}

    template <typename... Args>
    typename std::
        enable_if<!std::is_same<typename std::result_of<Fn(Args...)>::type,
                                void>::value,
                  typename std::result_of<Fn(Args...)>::type>::type
        apply(Args &&... args) {

        auto t0 = std::chrono::high_resolution_clock::now();
        auto result = fn_(std::forward<Args>(args)...);
        auto t1 = std::chrono::high_resolution_clock::now();

        time_ = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t1 - t0).count();

        return result;
    }

    template <typename... Args>
    typename std::
        enable_if<std::is_same<typename std::result_of<Fn(Args...)>::type,
                               void>::value,
                  void>::type
        apply(Args &&... args) {

        auto t0 = std::chrono::high_resolution_clock::now();
        fn_(std::forward<Args>(args)...);
        auto t1 = std::chrono::high_resolution_clock::now();

        time_ = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t1 - t0).count();
    }

    double time() const { return time_; }

  private:
    Fn fn_;
    double time_ = 0.0;
};

template <typename Fn>
FunctionTimer<Fn> make_timer(Fn &&fn) {
    return FunctionTimer<Fn>{std::forward<Fn>(fn)};
}

} // namespace utility
} // namespace tcc

#endif // TCC_UTILITY_FUNCTIONTIMER_H
