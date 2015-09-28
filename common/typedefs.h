#pragma once 
#ifndef TCC_COMMON_TYPEDEFS_H
#define TCC_COMMON_TYPEDEFS_H

#include "namespaces.h"
#include "foward_declerations.h"
#include "../include/eigen.h"

#include <type_traits>

// TileArray typedefs
using SpPolicy = TA::SparsePolicy;
using DnPolicy = TA::DensePolicy;

using TRange1 = TA::TiledRange1;
using TRange = TA::TiledRange;

using SpShapeF = TA::SparseShape<float>;

template<unsigned int DIM, typename Policy>
using TAArray = TA::Array<double, DIM, TA::TensorD, Policy>;

template <unsigned int DIM, typename Tile, typename Policy>
using DArray = TA::Array<double, DIM, Tile, Policy>;

// MPQC typedefs
template<typename E>
using Epool=tints::EnginePool<E>;

template<typename T>
using Tile = tcc::tensor::Tile<T>;

// Eig Typedefs
using MatrixD = Eig::Matrix<double, Eig::Dynamic, Eig::Dynamic, Eig::RowMajor>;
using VectorD = Eig::VectorXd;

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
