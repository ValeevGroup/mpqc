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

  class MOBlock {

  public:
    typedef std::pair<std::size_t, std::size_t> data_type;
    typedef std::pair<std::vector<std::size_t>, std::vector<std::size_t>> result_type;

    static char occ_char[4];
    static char vir_char[4];

    MOBlock(const tcc::TRange1Engine &tre) : occ_range_(0ul,
                                                        tre.get_occ_blocks()),
                                             vir_range_(tre.get_occ_blocks(),
                                                        tre.get_all_blocks()) { }

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
        if (std::find(occ_char, occ_char+4, i) != occ_char+4) {
          result.first.push_back(occ_range_.first);
          result.second.push_back(occ_range_.second);
        }
        else if (std::find(vir_char, vir_char+4, i) != vir_char+4) {
          result.first.push_back(vir_range_.first);
          result.second.push_back(vir_range_.second);
        } else {
          throw std::runtime_error("wrong key index");
        }
      }
      return result;
    }

  private:
    data_type occ_range_;
    data_type vir_range_;
  };

  char MOBlock::occ_char[4] = {'i', 'j', 'k', 'l'};
  char MOBlock::vir_char[4] = {'a', 'b', 'c', 'd'};
}
#endif //TILECLUSTERCHEM_MO_BLOCK_H
