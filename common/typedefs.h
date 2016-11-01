#pragma once
#ifndef MPQC_COMMON_TYPEDEFS_H
#define MPQC_COMMON_TYPEDEFS_H

#include "mpqc/math/external/eigen/eigen.h"

#include <type_traits>
#include <vector>
#include <memory>
#include "forward_declarations.h"

/////////////////////////////////////////////////////////
// TiledArray typedefs
/////////////////////////////////////////////////////////
using SpPolicy = TA::SparsePolicy;
using DnPolicy = TA::DensePolicy;

using TRange1 = TA::TiledRange1;
using TRange = TA::TiledRange;

using SpShapeF = TA::SparseShape<float>;

template <unsigned int, typename Tile, typename Policy>
using DArray = TA::DistArray<Tile, Policy>;

/////////////////////////////////////////////////////////
// MPQC typedefs
/////////////////////////////////////////////////////////
template <typename E>
using Epool = mpqc::integrals::EnginePool<E>;

template <typename T>
using Tile = mpqc::tensor::Tile<T>;

template <typename T>
using DecompTensor = mpqc::tensor::DecomposedTensor<T>;

using DecompTensorD = DecompTensor<double>;

/////////////////////////////////////////////////////////
// Eigen typedefs
/////////////////////////////////////////////////////////
using MatrixD = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using VectorD = Eigen::VectorXd;
using Vec3D = Eigen::Vector3d;

/////////////////////////////////////////////////////////
// Libint typedefs
/////////////////////////////////////////////////////////
using Shell = libint2::Shell;
using ShellVec = std::vector<libint2::Shell>;

/////////////////////////////////////////////////////////
// typetraits typedefs
/////////////////////////////////////////////////////////
template <typename T>
using remove_ref_t = typename std::remove_reference<T>::type;

template <typename T>
using remove_const_t = typename std::remove_const<T>::type;

template <typename T>
using remove_ref_const_t = typename std::remove_const<remove_ref_t<T>>::type;

template <typename T>
using result_of_t = typename std::result_of<T>::type;

template <bool B, typename U>
using enable_if_t = typename std::enable_if<B, U>::type;

template <typename T, typename ...Args>
std::unique_ptr<T> make_unique(Args &&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


#endif // MPQC_COMMON_TYPEDEFS_H
