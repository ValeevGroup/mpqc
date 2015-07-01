//
// Created by Chong Peng on 6/30/15.
//

#ifndef TILECLUSTERCHEM_TARRAY_BLOCK_H
#define TILECLUSTERCHEM_TARRAY_BLOCK_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

namespace tcc {

  template<typename T, unsigned int DIM,
          typename Tile, typename Policy, typename BlockEngine>
  class TArrayBlock {

  public:
    typedef TA::Array <T, DIM, Tile, Policy> TArray;

    /// Default constructor
    TArrayBlock() : array_(), block_engine_() { }

    /// Constructor
    TArrayBlock(const TArray &array,
                const std::shared_ptr<BlockEngine> block_engine)
            : array_(array), block_engine_(block_engine) { }

    /// Copy constructor

    /// shallow copy
    TArrayBlock(const TArrayBlock &other)
            : array_(other.array_), block_engine_(other.block_engine_) { }

    /// Assignment operator

    /// shallow copy
    TArrayBlock &operator=(const TArrayBlock &other) {
      array_ = other.array_;
      block_engine_ = other.block_engine_;

      return *this;
    }

    /// create expression

    /// \param vars A string with a comma-separated list of variables
    /// \return A non-const tensor block expression object
    TA::expressions::BlkTsrExpr<TArray>
    operator()(const std::string &vars) {

      std::pair<std::vector<std::size_t>, std::vector<std::size_t>> range =
              block_engine_->get(vars);
      return array_(vars).block(range.first, range.second);
    }

    /// Array accessor

    /// \return A reference to the array
    TArray &get_array() const {
      return &array_;
    }

  private:
    TArray array_;
    std::shared_ptr<BlockEngine> block_engine_;

  };

}
#endif //TILECLUSTERCHEM_TARRAY_BLOCK_H
