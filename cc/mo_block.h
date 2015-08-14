//
// Created by Chong Peng on 6/30/15.
//

#ifndef TILECLUSTERCHEM_MO_BLOCK_H
#define TILECLUSTERCHEM_MO_BLOCK_H

#include <memory>
#include <utility>
#include <vector>
#include "trange1_engine.h"

namespace tcc {

    // an easier interface to use TA's block expression to block by MO
    // use this class to initialize TArrayBlock class to block at occ or vir
    // occ: i,j,k,l
    // vir: a,b,c,d
    // all: p,q,r,s
    //
  class MOBlock {

  public:
    typedef std::pair<std::size_t, std::size_t> data_type;
    typedef std::pair<std::vector<std::size_t>, std::vector<std::size_t>> result_type;

    static const char occ_char[2];
    static const char vir_char[2];
    static const char all_char[2];

    /// Default constructor
    MOBlock() : occ_range_(), vir_range_(), all_range_() { }

    /// constructor

    MOBlock(const tcc::TRange1Engine &tre) : occ_range_(0ul, tre.get_occ_blocks()),
                                             vir_range_(tre.get_occ_blocks(),
                                                        tre.get_all_blocks()),
                                             all_range_(0ul, tre.get_all_blocks())
    { }

    /// Copy constructor

    /// deep copy
    MOBlock(const MOBlock &other) : occ_range_(other.occ_range_),
                                    vir_range_(other.vir_range_),
                                    all_range_(other.all_range_) { }

    /// Assignment operator

    /// deep copy
    MOBlock &operator=(const MOBlock &other) {
      occ_range_ = other.occ_range_;
      vir_range_ = other.vir_range_;
      all_range_ = other.all_range_;

      return *this;
    }

    /// get low_bound and upper_bound for TiledArray blocking

    /// vars is a expression string for TiledArray expression
    result_type get(const std::string &vars) {
      //select keys from vars
      std::vector<char> keys;
      for (std::size_t i = 0ul; i < vars.length(); ++i) {
        if (std::isalpha(vars[i])) {
          keys.push_back(vars[i]);
        }
      }

      //create ranges based on keys
      result_type result;
      for (auto i : keys) {
        if (i >= occ_char[0] && i <= occ_char[1]) {
          result.first.push_back(occ_range_.first);
          result.second.push_back(occ_range_.second);
        }
        else if (i >= vir_char[0] && i <= vir_char[1]) {
          result.first.push_back(vir_range_.first);
          result.second.push_back(vir_range_.second);
        }
        else if (i >= all_char[0] && i <= all_char[1]) {
          result.first.push_back(all_range_.first);
          result.second.push_back(all_range_.second);
        } else {
          throw std::runtime_error("wrong key index");
        }
      }
      return result;
    }

  private:
    data_type occ_range_;
    data_type vir_range_;
    data_type all_range_;
  };

  // set the range of key index to use
  const char MOBlock::occ_char[2] = {'i', 'l'};
  const char MOBlock::vir_char[2] = {'a', 'd'};
  const char MOBlock::all_char[2] = {'p', 's'};
}
#endif //TILECLUSTERCHEM_MO_BLOCK_H
