//
// string.h
//
// Copyright (C) 2012 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to_
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

#ifndef _mpqc_src_lib_chemistry_qc_nbody_string_h
#define _mpqc_src_lib_chemistry_qc_nbody_string_h

#include <bitset>
#include <vector>
#include <set>
#include <ostream>
#include <iostream>
#include <array>
#include <boost/dynamic_bitset.hpp>

namespace sc {

  /**
   * a "dense" string represents occupancies of a set of Ns states by a fixed-width bitstring
   * @tparam Ns the number of states
   */
  template <size_t Ns=64ul>
  class FermionOccupationNBitString : public std::bitset<Ns> {
    public:
      typedef size_t state_index_t;

      /**
       * Constructs an empty set of states
       */
      FermionOccupationNBitString() {
      }

      /**
       * Constructs FermionOccupationNBitString using a (possibly-empty) set of indices of occupied states
       * @param occupied_states
       */
      explicit FermionOccupationNBitString(const std::vector<state_index_t>& occupied_states) {
        for(auto v : occupied_states) {
          (*this)[v].flip();
        }
      }
      ~FermionOccupationNBitString() {}

      /**
       * are all states empty?
       * @return true if all states are empty
       */
      bool empty() const {
        return this->none();
      }

      /**
       * empties all states
       */
      void reset() {
        std::bitset<Ns>::reset();
      }

      /**
       * Reports the number of occupied states
       * @return the number of occupied states
       */
      size_t count() const {
        return std::bitset<Ns>::count();
      }

      /**
       * Removes a particle from_ state \a from_. Unsafe, but fast.
       * @param from_ the state from_ which the particle will be removed. The current status of the state is not checked.
       * @return reference to_ this string, for chaining operations (e.g. a.remove(0).remove(7).add(13) etc.)
       */
      FermionOccupationNBitString& remove(size_t from) {
        (*this)[from] = false;
        return *this;
      }

      /**
       * Adds a particle to_ state \a to_. Unsafe, but fast.
       * @param to_ the state to_ which the particle will be added. The current status of the state is not checked.
       * @return reference to_ this string, for chaining operations (e.g. a.remove(0).remove(7).add(13) etc.)
       */
      FermionOccupationNBitString& add(size_t to) {
        (*this)[to] = true;
        return *this;
      }

      /**
       * counts the number of bits in (pos1,pos2) (assumes pos2 >= pos1)
       * @param pos1
       * @param pos2
       * @return
       */
      size_t count(size_t pos1, size_t pos2) const {
        // to_ determine the number of particles crossed count the number of particles between to_ and from_
        FermionOccupationNBitString tmp(*this);
        tmp >>= (pos1+1);
        tmp <<= (Ns - (pos2 - pos1 - 1));
        return tmp.count();
      }

    private:
      static const size_t nstates_ = Ns;
  };

  template <class CharT, class Traits, size_t Ns>
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os,
             const FermionOccupationNBitString<Ns>& x)
  {
      os << std::bitset<Ns>(x);
      return os;
  }

  /**
   * a "dense" string represents occupancies of a set of Ns states by a bitstring
   * @tparam Ns the number of states
   */
  class FermionOccupationDBitString : public boost::dynamic_bitset<> {
    public:
      typedef size_t state_index_t;

      /**
       * Constructs an empty set of N states
       */
      FermionOccupationDBitString(size_t N) : boost::dynamic_bitset<>(N, 0ul), nstates_(N) {
      }

      /**
       * Constructs FermionOccupationDBitString using a (possibly-empty) set of indices of occupied states
       * @param occupied_states
       */
      explicit FermionOccupationDBitString(size_t N, const std::vector<state_index_t>& occupied_states) :
        boost::dynamic_bitset<>(N, 0ul), nstates_(N)
      {
        for(auto v : occupied_states) {
          (*this)[v].flip();
        }
      }
      ~FermionOccupationDBitString() {}

      /**
       * are all states empty?
       * @return true if all states are empty
       */
      bool empty() const {
        return this->none();
      }

      /**
       * empties all states
       */
      void reset() {
        boost::dynamic_bitset<>::reset();
      }

      /**
       * Reports the number of occupied states
       * @return the number of occupied states
       */
      size_t count() const {
        return boost::dynamic_bitset<>::count();
      }

      /**
       * Removes a particle from_ state \a from_. Unsafe, but fast.
       * @param from_ the state from_ which the particle will be removed. The current status of the state is not checked.
       * @return reference to_ this string, for chaining operations (e.g. a.remove(0).remove(7).add(13) etc.)
       */
      FermionOccupationDBitString& remove(size_t from) {
        (*this)[from] = false;
        return *this;
      }

      /**
       * Adds a particle to_ state \a to_. Unsafe, but fast.
       * @param to_ the state to_ which the particle will be added. The current status of the state is not checked.
       * @return reference to_ this string, for chaining operations (e.g. a.remove(0).remove(7).add(13) etc.)
       */
      FermionOccupationDBitString& add(size_t to) {
        (*this)[to] = true;
        return *this;
      }

      /**
       * counts the number of bits in (pos1,pos2) (assumes pos2 >= pos1)
       * @param pos1
       * @param pos2
       * @return
       */
      size_t count(size_t pos1, size_t pos2) const {
        // to_ determine the number of particles crossed count the number of particles between to_ and from_
        FermionOccupationDBitString tmp(*this);
        tmp >>= (pos1+1);
        tmp <<= (nstates_ - (pos2 - pos1 - 1));
        return tmp.count();
      }

    private:
      size_t nstates_;

  };

  template <class CharT, class Traits>
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os,
             const FermionOccupationDBitString& x)
  {
      os << boost::dynamic_bitset<>(x);
      return os;
  }

#if 1
  /**
   * a block-"sparse" string represents occupancies of an arbitrarily-large set of states as a set of alternating unoccupied/occupied blocks
   */
  class FermionOccupationBlockString {
    public:
      typedef size_t state_index_t;

      /// represents a continuous block of states of same occupancy
      struct Block {
        Block(size_t o, bool n = true) : offset(o), length(1ul), occ(n) {}
        Block(size_t o, size_t l, bool n = true) : offset(o), length(l), occ(n) {}
        size_t offset;
        mutable size_t length;
        bool occ;
        bool operator<(const Block& other) const { return offset < other.offset; }
      };
      typedef std::set<Block> Blocks; //< stores blocks

      /**
       * Constructs an empty set of Ns states
       * @param Ns number of states
       */
      explicit FermionOccupationBlockString(size_t Ns) : nstates_(Ns), blocks_() {
        blocks_.insert(Block(0,Ns,false));
      }

      /**
       * Constructs FermionOccupationBlockString using a (possibly-empty) set of indices of occupied states
       * @param occupied_states
       */
      explicit FermionOccupationBlockString(size_t Ns, const std::vector<state_index_t>& occupied_states) : nstates_(Ns), blocks_() {
        blocks_.insert(Block(0,Ns,false));
        // not efficient -- should find blocks
        for(auto v : occupied_states) {
          this->add(v);
        }
      }

      /**
       * Constructs FermionOccupationBlockString using a (possibly-empty) set of indices of occupied states
       * @param occupied_states
       */
      explicit FermionOccupationBlockString(size_t Ns, Blocks& blocks) : nstates_(Ns), blocks_() {
        std::swap(blocks_,blocks);
      }

      ~FermionOccupationBlockString() {}

      /**
       * are all states empty?
       * @return true if all states are empty
       */
      bool empty() const {
        return this->count() == 0;
      }

      /**
       * empties all states
       */
      void reset() {
        blocks_.clear();
        blocks_.insert(Block(0,nstates_,false));
      }

      /**
       * Reports the total number of states
       * @return the total number states
       */
      size_t size() const {
        return nstates_;
      }

      /**
       * Reports the number of occupied states
       * @return the number of occupied states
       */
      size_t count() const {
        size_t c = 0;
        for(auto v : blocks_) {
          if (v.occ) c += v.length;
        }
        return c;
      }

      /**
       * Returns the occupancy of state i
       * @param i the state
       * @return true is the state is occupied
       */
      bool operator[](size_t i) const {
        auto ub_iter = blocks_.upper_bound(Block(i));
        auto result_iter = --ub_iter;
        return result_iter->occ;
      }

      /**
       * Removes a particle from_ state \a from_. Unsafe, but fast.
       * @param from_ the state from_ which the particle will be removed. The current status of the state is not checked.
       * @return reference to_ this string, for chaining operations (e.g. a.remove(0).remove(7).add(13) etc.)
       */
      FermionOccupationBlockString& remove(size_t from) {

        // 4 possibilities:
        // 1) state is at the beginning/end of a block -- shrink the occupied block and extend an unoccupied block
        // 2) state is in the middle of a block -- split the occupied block, create an unoccupied block
        // 3) state is the only one in the block -- delete the occupied block, merge the unoccupied blocks
        auto ub_iter = blocks_.upper_bound(Block(from));
        auto result_iter = --ub_iter;
        assert(result_iter->occ == true);

        if (result_iter->length == 1) { // only
          auto uocc_block = this->block_erase(result_iter);
        }
        else if (from == result_iter->offset) { // front
          auto occ_block = this->block_pop_front(result_iter);
          if (occ_block != blocks_.begin()) {
            auto uocc_block = this->block_push_back(--occ_block);
          }
          else
            auto uocc_block = blocks_.insert(occ_block, Block(0, 1, false));
        }
        else if (from == result_iter->offset + result_iter->length - 1) { // back
          auto occ_block = this->block_pop_back(result_iter);
          if (++occ_block != blocks_.end()) {
            auto uocc_block = this->block_push_front(occ_block);
          }
          else
            auto uocc_block = blocks_.insert(occ_block, Block(0, 1, false));
        }
        else { // in the middle
          this->block_split(result_iter, from);
        }

        return *this;
      }

      /**
       * Adds a particle to_ state \a to_. Unsafe, but fast.
       * @param to_ the state to_ which the particle will be added. The current status of the state is not checked.
       * @return reference to_ this string, for chaining operations (e.g. a.remove(0).remove(7).add(13) etc.)
       */
      FermionOccupationBlockString& add(size_t to) {

        // 4 possibilities:
        // 1) state is at the beginning/end of a block -- shrink the occupied block and extend an unoccupied block
        // 2) state is in the middle of a block -- split the occupied block, create an unoccupied block
        // 3) state is the only one in the block -- delete the occupied block, merge the unoccupied blocks
        auto ub_iter = blocks_.upper_bound(Block(to));
        auto result_iter = --ub_iter;
        assert(result_iter->occ == false);

        if (result_iter->length == 1) { // only
          auto occ_block = this->block_erase(result_iter);
        }
        else if (to == result_iter->offset) { // front
          auto uocc_block = this->block_pop_front(result_iter);
          if (uocc_block != blocks_.begin()) {
            auto occ_block = this->block_push_back(--uocc_block);
          }
          else
            auto occ_block = blocks_.insert(uocc_block, Block(0, 1, true));
        }
        else if (to == result_iter->offset + result_iter->length - 1) { // back
          auto uocc_block = this->block_pop_back(result_iter);
          if (++uocc_block != blocks_.end()) {
            auto occ_block = this->block_push_front(uocc_block);
          }
          else
            auto occ_block = blocks_.insert(uocc_block, Block(nstates_-1, 1, true));
        }
        else { // in the middle
          this->block_split(result_iter, to);
        }

        return *this;
      }

      /**
       * counts the number of bits in (pos1,pos2) (assumes pos2 >= pos1)
       * @param pos1
       * @param pos2
       * @return
       */
      size_t count(size_t pos1, size_t pos2) const {
        // to_ determine the number of particles crossed count the number of particles between to_ and from_
        assert(pos1 <= pos2);
        auto ub1_iter = blocks_.upper_bound(Block(pos1));
        auto blk1_iter = --ub1_iter;
        auto ub2_iter = blocks_.upper_bound(Block(pos2));
        auto blk2_iter = --ub2_iter;
        size_t c = (blk1_iter->occ ? (blk1_iter->offset + blk1_iter->length - pos1 - 1) : 0);
        c += (blk2_iter->occ ? (pos2 - blk2_iter->offset) : 0);
        ++blk1_iter; --blk2_iter;
        for(auto i=blk1_iter; i!=blk2_iter; ++i)
          if (i->occ)
            c += i->length;
        std::cout << "count = " << c << std::endl;
        return c;
      }

      boost::dynamic_bitset<> to_bitset() const {
        boost::dynamic_bitset<> result(nstates_);
        for(auto v : blocks_) {
          if (v.occ)
            for(size_t l=0, pos=v.offset;
                l<v.length; ++l, ++pos) {
              result.flip(pos);
            }
        }
        return result;
      }

      /// XORs two strings
      FermionOccupationBlockString operator^(const FermionOccupationBlockString& other) const {
        assert(size() == other.size());
        Blocks result_blocks;
        auto w = result_blocks.begin();
        auto u = other.blocks().begin();
        auto v = blocks().begin();
        auto red = v;
        auto black = u;
        size_t red_end = red->offset+red->length;
        size_t black_end = black->offset+black->length;
        const size_t size2 = size() * 2;
        while (red_end + black_end != size2) {
          if (red_end <= black_end) {
            w = result_blocks.insert(w, Block(red->offset, red->length, red->occ ^ black->occ));
            ++red;
            red_end = red->offset+red->length;
            if (red->offset == black_end) {
              ++black;
              black_end = black->offset+black->length;
            }
          }
          else {
            w = result_blocks.insert(w, Block(red->offset, black_end - red->offset, red->occ ^ black->occ));
            ++black;
            black_end = black->offset+black->length;
            std::swap(red,black);
            std::swap(red_end,black_end);
          }
        }
        w = result_blocks.insert(w, Block(red->offset, red->length, red->occ ^ black->occ));
        FermionOccupationBlockString result(size(), result_blocks);
        return result;
      }

    private:
      size_t nstates_;
      Blocks blocks_;
      Blocks& blocks() { return blocks_; }
      const Blocks& blocks() const { return blocks_; }

      Blocks::const_iterator block_pop_front(Blocks::const_iterator block) {
        Block popped_block(*block);
        --popped_block.length;
        ++popped_block.offset;
        blocks_.erase(block);
        const auto result_iter = blocks_.insert(popped_block).first;
        return result_iter;
      }
      Blocks::const_iterator block_pop_back(Blocks::const_iterator block) {
        --(block->length);
        return block;
      }
      Blocks::const_iterator block_push_front(Blocks::const_iterator block) {
        Block popped_block(*block);
        ++popped_block.length;
        --popped_block.offset;
        blocks_.erase(block);
        const auto result_iter = blocks_.insert(popped_block).first;
        return result_iter;
      }
      Blocks::const_iterator block_push_back(Blocks::const_iterator block) {
        ++(block->length);
        return block;
      }

      Blocks::const_iterator block_erase(Blocks::const_iterator block) {
        const bool have_prev = (block != blocks_.end());
        Blocks::const_iterator next_block(block); ++next_block;
        const bool have_next = (next_block != blocks_.end());

        if (have_prev) {
          Blocks::const_iterator prev_block(block); --prev_block;
          prev_block->length += block->length;
          blocks_.erase(block);
          if (have_next) {
            prev_block->length += next_block->length;
            blocks_.erase(next_block);
          }
          return prev_block;
        }
        else if (have_next) { // not have_prev
          Block new_block(*next_block);
          const size_t length = block->length;
          new_block.offset -= length;
          new_block.length += length;
          blocks_.erase(block);
          blocks_.erase(next_block);
          auto result = blocks_.insert(new_block).first;
          return result;
        }
        else { // erasing the only block
          blocks_.clear();
          blocks_.insert(Block(0,nstates_,false));
          return blocks_.begin();
        }
        assert(false); // unreachable
        return blocks_.begin();
      }

      Blocks::const_iterator block_split(Blocks::const_iterator block,
                                                  size_t pos) {
        const size_t orig_length = block->length;
        const size_t offset = block->offset;
        block->length = (pos - block->offset);

        Block split_block(pos, 1, !block->occ);
        auto iter = blocks_.insert(block, split_block);

        Block remainder_block(pos+1, offset+orig_length-pos-1, block->occ);
        blocks_.insert(iter, remainder_block);

        return iter;
      }

  };

  template <class CharT, class Traits>
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os,
             const FermionOccupationBlockString& x)
  {
      os << x.to_bitset();
      return os;
  }

#endif

  /**
   * basic Nb-body number-conserving (nc) operator in sp representation
   * @tparam Nb rank of the operator (number of particles that possibly change states)
   */
  template <size_t Nb, typename FString>
  class FermionBasicNCOper {
    public:
      template <typename Int> FermionBasicNCOper(Int t, Int f) : to_(t) , from_(f) {
      }

      /**
       * reports the state to which this operator places a particle
       * @return the state to which this operator places a particle
       */
      size_t to() const { return to_; }

      /**
       * reports the state from which this operator removes a particle
       * @return the state from which this operator removes a particle
       */
      size_t from() const { return from_; }

      /**
       * applies this operator to_ FermionOccupationNBitString os. to_==from_ is allowed.
       * @param os the FermionOccupationNBitString object
       * @return FermionOccupationNBitString obtained by removing a particle from_ \a from_ and adding a particle to_ \a to_.
       */
      void apply(FString& os) const {
        if (os[from_] == true && (os[to_] == false || from_ == to_)) {
          os.remove(from_).add(to_);
        }
        else
          os.reset();
      }

      /**
       * same as apply(), but returns whether application operator changes the sign of the state; the sign changes if
       * the number of occupied states "crossed" by the operator is odd
       * @param os the FermionOccupationNBitString object
       * @return true is the sign of the state changes
       */
      bool apply_sign(FString& os) const {
        if (os[from_] == true && (os[to_] == false || from_ == to_)) {
          os.remove(from_).add(to_);
          // to_ determine the number of particles crossed, count the number of particles between to_ and from_
          if (to_ > from_) {
            const bool sign_changes = os.count(from_,to_) & size_t(1);
            return sign_changes;
          } else { // to_ <= from_
            const bool sign_changes = os.count(to_,from_) & size_t(1);
            return sign_changes;
          }
        }
        else {
          os.reset();
          return false;
        }
      }

      /**
       * similar to_ apply_sign(), but keeps the argument unchanged
       * @param os the FermionOccupationNBitString object
       * @return new pair[sign_change,FermionOccupationNBitString]
       */
      std::pair<bool,FString > operator()(const FString& os) const {
        FString result(os);
        const bool sign_changes = this->apply_sign(result);
        return std::make_pair(sign_changes,result);
      }

    private:
      size_t to_, from_;
  };

  template <class CharT, class Traits, size_t Nb, typename FString>
  std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os,
             const FermionBasicNCOper<Nb, FString>& o)
  {
      os << o.to() << "<-" << o.from();
      return os;
  }

} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
