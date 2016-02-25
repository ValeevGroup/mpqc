#pragma once
#ifndef MPQC_UTILITY_PRINTSIZEINFO_H
#define MPQC_UTILITY_PRINTSIZEINFO_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"

#include "array_info.h"

#include <string>

namespace mpqc {
namespace utility {

template<typename Array>
void print_size_info(Array const &A, std::string const &name){
   auto &world = A.get_world(); 
   auto sizes = array_storage(A);

   if(world.rank() == 0){
       std::cout << "Printing size information for " << name << "\n";


       std::cout << "\tFull     = " << sizes[0] << " GB\n" 
                 << "\tSparse   = " << sizes[1] << " GB\n"
                 << "\tLow Rank = " << sizes[2] << " GB\n\n";
   }

   world.gop.fence();
}

template<typename Array>
void wprint_size_info(Array const &A, std::wstring const &name){
    auto &world = A.get_world();
    auto sizes = array_storage(A);

    if(world.rank() == 0){
        std::cout << "Printing size information for ";
        std::wcout << name << "\n";


        std::cout << "\tFull     = " << sizes[0] << " GB\n"
        << "\tSparse   = " << sizes[1] << " GB\n"
        << "\tLow Rank = " << sizes[2] << " GB\n\n";
    }

    world.gop.fence();
}
} // namespace utility
} // namespace mpqc

#endif // MPQC_UTILITY_PRINTSIZEINFO_H
