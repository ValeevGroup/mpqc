#pragma once
#ifndef TCC_PURIFICATION_PURIFICATIONDEVEL_H
#define TCC_PURIFICATION_PURIFICATIONDEVEL_H

#include <cstdlib>

namespace tcc {
namespace pure {

template <typename Array>
double compute_trace(Array const &array) {
    // calculate the trace of the array
    double trace = 0;
    return trace;
}

template <typename Array>
Array initial_guess(Array const &array) {
    Array guess(array.get_world(), array.trange());
    guess.set_all_local(1.0);
    return guess;
}

template <typename Array>
void trace_resetting_poly(Array &array, std::size_t occ, double trace) {
    array("i,j") = array("i,j") * array("j,i");
    array.get_world().gop.fence();
}


class purifier {
  public:
    explicit purifier(double cut = 1e-06) : cut_{cut} {}

    template <typename Array, typename Polynomial> 
    void operator()(Array const &array, std::size_t occupation,
                    Polynomial poly = trace_resetting_poly<Array>){
        auto P = initial_guess(array);
        auto trace = compute_trace(P);
        poly(P, occupation, trace);
    }

  private:
    double cut_;
}; // class purifier

} // namespace pure
} // namespace tcc

#endif /* end of include guard: TCC_PURIFICATION_PURIFICATIONDEVEL_H */
