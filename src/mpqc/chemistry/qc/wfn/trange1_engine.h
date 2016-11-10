//
// Created by Chong Peng on 6/25/15.
//

#ifndef MPQC_BLOCK_SIZE_ENGINE_H
#define MPQC_BLOCK_SIZE_ENGINE_H

#include <tiledarray.h>

#include <algorithm>

namespace mpqc {
class TRange1Engine {
 public:
  TRange1Engine()
      : occ_(0ul),
        all_(0ul),
        vir_(0ul),
        occ_block_size_(0ul),
        vir_block_size_(0ul),
        nfrozen_(0ul),
        tr_active_occupied_(),
        tr_virtual_(),
        tr_all_() {}

  TRange1Engine(const std::size_t occ, const std::size_t all,
                const std::size_t occ_block_size,
                const std::size_t vir_block_size,
                const std::size_t nfrozen = 0ul)
      : occ_(occ),
        all_(all),
        vir_(all - occ),
        occ_block_size_(std::min(occ_block_size, occ)),
        vir_block_size_(std::min(vir_block_size, vir_)),
        nfrozen_(nfrozen) {
    init();
  }

  void init();

  TA::TiledRange1 compute_range(std::size_t range, std::size_t block_size);

  /// get TiledRange1 for active occuppied space
  TA::TiledRange1 get_active_occ_tr1() const { return tr_active_occupied_; }

  /// get TiledRange1 for occuppied space
  TA::TiledRange1 get_occ_tr1() const { return tr_occupied_; }

  /// get TiledRange1 for active virtual occuppied space
  TA::TiledRange1 get_vir_tr1() const { return tr_virtual_; }

  /// /// get TiledRange1 for all obs space
  TA::TiledRange1 get_all_tr1() const { return tr_all_; }

  std::size_t get_occ() const { return occ_; }

  /// get the number of active occ orbitals
  std::size_t get_active_occ() const { return occ_ - nfrozen_; }

  std::size_t get_vir() const { return vir_; }

  std::size_t get_all() const { return all_; }

  /// get the number of active orbitals
  std::size_t get_active_all() const { return all_ - nfrozen_; }

  std::size_t get_nfrozen() const { return nfrozen_; }

  std::size_t get_occ_block_size() const { return occ_block_size_; }

  std::size_t get_vir_block_size() const { return vir_block_size_; }

  /// get the number of blocks in active occupied
  std::size_t get_active_occ_blocks() const { return active_occ_blocks_; }

  /// get the number of blocks in virtual
  std::size_t get_vir_blocks() const { return vir_blocks_; }

 private:
  TA::TiledRange1 tr_occupied();

  TA::TiledRange1 tr_active_occupied();

  TA::TiledRange1 tr_virtual();

  TA::TiledRange1 tr_all();

 private:
  /// number of occupied orbitals
  std::size_t occ_;

  /// number of all orbitals
  std::size_t all_;

  /// number of virtual orbitals
  std::size_t vir_;

  std::size_t occ_block_size_;
  std::size_t vir_block_size_;

  /// number of frozen orbitals
  std::size_t nfrozen_;

  TA::TiledRange1 tr_active_occupied_;
  TA::TiledRange1 tr_occupied_;
  TA::TiledRange1 tr_virtual_;
  TA::TiledRange1 tr_all_;

  /// number of active occupied blocks
  std::size_t active_occ_blocks_;

  /// number of virtual blocks
  std::size_t vir_blocks_;
};
}

#endif  // MPQC_BLOCK_SIZE_ENGINE_H
