/*
 *  This file is a part of TiledArray.
 *  Copyright (C) 2017  Virginia Tech
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Justus Calvin
 *  Department of Chemistry, Virginia Tech
 *
 *  dist_array.cpp
 *  Feb 5, 2016
 *
 */

#include <TiledArray/dist_eval/array_eval.h>
#include <TiledArray/dist_eval/contraction_eval.h>
#include <TiledArray/tile_op/contract_reduce.h>
#include <TiledArray/tile_op/noop.h>
#include <TiledArray/tile_op/shift.h>
#include <TiledArray/tile_op/unary_wrapper.h>

namespace TiledArray {
namespace detail {

template class Summa<
    DistEval<LazyArrayTile<TiledArray::Tensor<double>,
                           UnaryWrapper<TiledArray::Noop<
                               TiledArray::Tensor<double>, true> > >,
             TiledArray::SparsePolicy>,
    DistEval<LazyArrayTile<TiledArray::Tensor<double>,
                           UnaryWrapper<TiledArray::Shift<
                               TiledArray::Tensor<double>, false> > >,
             TiledArray::SparsePolicy>,
    TiledArray::ContractReduce<TiledArray::Tensor<double>,
                               TiledArray::Tensor<double>, double>,
    TiledArray::SparsePolicy>;

}  // namespace detail
}  // namespace TiledArray
