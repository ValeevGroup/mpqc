#include "../include/btas.h"
#include "../include/tiledarray.h"
#include "gtest.h"

#include <iostream>

TEST(BTAS, DefaultConstructor){
    btas::Tensor<double> a(2,2,2);
    EXPECT_EQ(3, a.rank()); 
    EXPECT_EQ(8, a.size());
    
    a.fill(1.0);
    auto check_array = 
        [&](btas::Tensor<double> const &array) {
            for(auto const& e : array){
                if(e != 1.0) {
                    return false;
                }
            }
            return true;
        };
    EXPECT_TRUE(check_array(a)) << "At least one element was not 1.0.";
}

