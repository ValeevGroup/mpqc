
#include <TiledArray/contraction_eval.h>

#ifndef TILEDARRAY_HEADER_ONLY
extern template class Summa<DistEval<LazyArrayTile<TiledArray::Tensor<double>, UnaryWrapper<TiledArray::Noop<TiledArray::Tensor<double>, true> > >, TiledArray::SparsePolicy>, DistEval<LazyArrayTile<TiledArray::Tensor<double>, UnaryWrapper<TiledArray::Shift<TiledArray::Tensor<double>, false> > >, TiledArray::SparsePolicy>, TiledArray::ContractReduce<TiledArray::Tensor<double>, TiledArray::Tensor<double>, double>, TiledArray::SparsePolicy>;
#endif  // !defined(TILEDARRAY_HEADER_ONLY)
