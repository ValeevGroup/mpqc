#include "common/namespaces.h"
#include "common/typedefs.h"
#include "include/tiledarray.h"

#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/external/rapidjson/rapidjson.h>
#include <cereal/archives/json.hpp>

#include <utility>
#include <iostream>
#include <sstream>

int main(int argc, char *argv[]) {
    std::stringstream ss;
    {
        std::vector<double> my_vec = {1.0, 2.0, 3.0};
        std::array<double,4> my_arr = {{4.0, 5.0, 6.0, 7.0}};
        cereal::JSONOutputArchive oarchive_cout(std::cout);
        cereal::JSONOutputArchive oarchive(ss);
        oarchive_cout(CEREAL_NVP(my_vec));
        oarchive_cout(CEREAL_NVP(my_arr));
        std::cout << std::endl;
        oarchive(my_vec);
        oarchive(my_arr);
    }
    {
        cereal::JSONInputArchive iarchive(ss);
        std::vector<double> my_vec;
        std::array<double,4> my_arr;

        iarchive(my_vec, my_arr);

        std::cout << "\nVector = ";
        for(auto elem : my_vec){
            std::cout << elem << " ";
        }
        std::cout << "\nArray = ";
        for(auto elem : my_arr){
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
