#pragma once 
#ifndef TCC_COMMON_TYPEDEFS_H
#define TCC_COMMON_TYPEDEFS_H

#include "namespaces.h"
#include "foward_declerations.h"

#include <type_traits>

// Typedefs for my most commonly used arrays
template<unsigned int DIM, typename Policy>
using LRArray = TA::Array<double, DIM, tcc::tensor::TilePimpl<double>, Policy>;

template<unsigned int DIM, typename Policy>
using TAArray = TA::Array<double, DIM, TA::TensorD, Policy>;

template<typename T>
using Tile = tcc::tensor::Tile<T>;

// Useful typedefs for removing qualifiers from types.
template<typename T>
using remove_ref_t = typename std::remove_reference<T>::type; 

template<typename T>
using remove_const_t = typename std::remove_const<T>::type; 

template<typename T>
using remove_ref_const_t = typename std::remove_const<remove_ref_t<T>>::type; 

template<typename T>
using result_of_t = typename std::result_of<T>::type;

template<bool B, typename U>
using enable_if_t = typename std::enable_if<B,U>::type;

#endif // TCC_COMMON_TYPEDEFS_H
