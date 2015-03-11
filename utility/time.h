#pragma once
#ifndef TCC_UTILITY_TIME_H
#define TCC_UTILITY_TIME_H

#include "../common/typedefs.h"

#include <chrono>

namespace tcc {
namespace utility {
namespace time {

using t_point = std::chrono::high_resolution_clock::time_point;

t_point now() { return std::chrono::high_resolution_clock::now(); }

double duration_in_s(t_point const &t0, t_point const &t1) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0)
        .count();
}


// Class to time a function or function object.  Call apply to run and time
// the function.
template <typename Fn>
class FunctionTimer {
  public:
    FunctionTimer(Fn fn) : fn_{fn} {}

    template <typename... Args>
    enable_if_t<!std::is_same<result_of_t<Fn(Args...)>, void>::value,
                result_of_t<Fn(Args...)>>
    apply(Args &&... args)  {

        auto t0 = now();
        auto result = fn_(std::forward<Args>(args)...);
        auto t1 = now();

        time_ = duration_in_s(t0, t1);

        return result;
    }

    template <typename... Args>
    enable_if_t<std::is_same<result_of_t<Fn(Args...)>, void>::value, void>
    apply(Args &&... args) {

        auto t0 = now();
        fn_(std::forward<Args>(args)...);
        auto t1 = now();

        time_ = duration_in_s(t0, t1);
    }

    double time() const { return time_; }

  private:
    Fn fn_;
    double time_ = 0.0;
}; // FunctionTimer

template <typename Fn>
time::FunctionTimer<Fn> make_timer(Fn &&fn) {
    return time::FunctionTimer<Fn>{std::forward<Fn>(fn)};
}

} // namespace time


} // namespace utility

namespace tcc_time = utility::time;

} // namespace tcc

#endif // TCC_UTILITY_TIME_H
