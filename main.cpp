#include "utility/expand_container.h"

#include <iostream>
#include <array>

double func(std::string str, double x, double y) {
    std::cout << str << " ";
    return x + y;
}

double
func(std::string str, double x, double y, double z) {
    std::cout << str << " ";
    return x + y + z;
}

struct call_func {
    using return_type = double;
    template <typename ...Args>
    double operator()(std::string str, Args ...args){
        return func(str, args...);
    }
};

int main(int argc, char *argv[]) {
    std::array<double, 2> xy = {{1.0, 2.0}};
    std::string op = "The sum of x + y = ";
    auto caller_2 = tcc::utility::make_tuple_expander(xy, op);
    auto sum = caller_2(call_func{});
    std::cout << sum << std::endl;

    std::array<double, 3> xyz = {{1.0, 2.0, 3.0}};
    op = "The sum of x + y + z = ";
    auto caller_3 = tcc::utility::make_tuple_expander(xyz, op);
    sum = caller_3(call_func{});
    std::cout << sum << std::endl;
    return 0;
}
