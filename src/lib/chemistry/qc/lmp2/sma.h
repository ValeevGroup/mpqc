
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_sma_h
#define _chemistry_qc_lmp2_sma_h

//#define USE_HASH
//#define USE_STL_MULTIMAP 1

#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/integral.h>
#include <util/misc/autovec.h>
#include <util/misc/scint.h>
#include <util/misc/exenv.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/misc/regtime.h>
#include <util/misc/scexception.h>

#ifndef TWO_INDEX_SPECIALIZATIONS
#  define THREE_INDEX_SPECIALIZATIONS 1
#endif

#ifndef THREE_INDEX_SPECIALIZATIONS
#  define THREE_INDEX_SPECIALIZATIONS 1
#endif

#ifndef FOUR_INDEX_SPECIALIZATIONS
#  define FOUR_INDEX_SPECIALIZATIONS 1
#endif

#ifndef USE_STL_MULTIMAP
#  define USE_STL_MULTIMAP 0
#endif

#if ! USE_STL_MULTIMAP
#  include <chemistry/qc/lmp2/avlmmap.h>
#endif

#include <math.h>
#include <float.h>
#include <vector>
#include <map>
#ifdef USE_HASH
#  include <ext/hash_map>
#endif
#include <set>
#include <stdexcept>
#include <string>
#ifdef WORDS_BIGENDIAN
#  define r2(i) (i)
#  define r4(i) (i)
#else
// For little endian machines, the order of indices in BlockInfo<4>
// must be reversed so a fast comparison on an uint64_t can be used
// to order blocks.
// BlockInfo<3> uses r4 and zeros out the last byte.
// BlockInfo<2> uses r2.
//#  define r4(i) (3-(i))
#  define r2(i) ((~(i))&1)
#  define r4(i) ((~(i))&3)
#endif

namespace sc {

namespace sma2 {

    /** \brief An Range represent a set of integers, [0, N).
        For example, one Range could represent the atomic orbitals and
        another could represent the molecular orbitals.  An index is
        divided into user defined blocks. */
    class Range {
      private:
        int nindex_;
        std::vector<int> block_size_;
        std::vector<int> block_offset_;
        std::vector<int> index_to_block_;
        std::vector<int> function_order_;
        std::vector<int> range_order_;

        int max_block_size_;

        void init_offsets();
        void init_extent_blocking(const sc::Ref<sc::GaussianBasisSet> &bs);
        void init_order(int nbasis);
      public:
        enum BlockingMethod { AtomBlocking,
                              ShellBlocking,
                              FunctionBlocking,
                              ExtentBlocking };

        Range(const sc::Ref<sc::GaussianBasisSet> &,
              BlockingMethod b = ShellBlocking, int blocksize = 1);
        Range(int nindex, int block_size);
        Range(const std::vector<int> & block_size);
        Range();

        void init(const sc::Ref<sc::GaussianBasisSet> &,
                  BlockingMethod b = ShellBlocking, int blocksize = 1);
        void init(int nindex, int block_size);
        void clear();

        /// Returns the dimension of this range.
        int nindex() const { return nindex_; }
        /// Returns the number of blocks.
        int nblock() const { return block_size_.size(); }
        /// Returns the size for the given block.
        int block_size(int i) const { return block_size_[i]; }
        /// Returns the maximum block size.
        int max_block_size() const { return max_block_size_; }
        /// Returns the offset for the given block.
        int block_offset(int i) const { return block_offset_[i]; }
        /// Given an index, return the block in which that index is found.
        int index_to_block(int i) const { return index_to_block_[i]; }
        /** Given an index, return its offset within the block in which
            that index is found. */
        int index_to_offset(int i) const {
          return i - block_offset_[index_to_block_[i]];
        }
        /** Maps a basis function number to the number within this
            Range object.  This arises since ExtentBlocking must
            reorder shells to group shells with similar extents in
            the same block. */
        int basis_index_to_range_index(int i) const;
        int range_index_to_basis_index(int i) const;
        bool operator == (const Range &r) const;
        bool operator != (const Range &r) const { return ! operator == (r); }
        /// Returns true if all blocks are size 1.
        bool all_size_one() const;
        void print(std::ostream&o=sc::ExEnv::out0()) const;
        void write(sc::StateOut&) const;
        void read(sc::StateIn&);
    };

    /** \brief An IndexList is a vector of indices.  Each element
        corresponds to the position of an index in an array.  For
        example, the index list [3, 2] corresponds to the third and
        second indices, in that order, in a multidimensional array. */
    class IndexList {
      private:
        std::vector<int> indices_;
      public:
        IndexList();
        IndexList(int);
        IndexList(int,int);
        IndexList(int,int,int);
        IndexList(int,int,int,int);
        IndexList(const IndexList &);
        /** Append the indices from both IndexList objects, il1 and il2, to
            form this IndexList. */
        IndexList(const IndexList &il1, const IndexList &il2);
        IndexList(const std::vector<int> &);
        IndexList reverse_mapping() const;
        void append_additional_indices(const IndexList &);
        int &i(int a) { return indices_[a]; }
        const int &i(int a) const { return indices_[a]; }
        int n() const { return indices_.size(); }
        void set_n(int n) { indices_.resize(n); }
        static IndexList identity(int n);
        /// Returns true if this is 0, 1, ..., n-1.
        bool is_identity() const;
        /// Returns true if this is a permutation of 0, 1, ..., n-1.
        bool is_identity_permutation() const;
        void print(std::ostream &o=sc::ExEnv::outn()) const;
        bool operator>(const IndexList&il)const{return indices_>il.indices_;}
        bool operator<(const IndexList&il)const{return indices_<il.indices_;}
        bool operator==(const IndexList&il)const{return indices_==il.indices_;}
        bool operator!=(const IndexList&il)const{return indices_!=il.indices_;}
        bool operator>=(const IndexList&il)const{return indices_>=il.indices_;}
        bool operator<=(const IndexList&il)const{return indices_<=il.indices_;}
    };
    std::ostream &operator << (std::ostream &, const IndexList &);

    typedef unsigned short bi_t; // Block index type.

    /** \brief BlockInfo stores info about a block of data.
        This info includes the block numbers for each index. */
    template <int N>
    class BlockInfo {
      private:
#ifdef USE_BOUND
        mutable double bound_;
#endif
        bi_t block_[N];
      public:
        BlockInfo() {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        /** Initialize the blocks to the values in the argument, v.  The
            number of blocks assigned is the smaller of N and v.size(). */
        BlockInfo(const std::vector<bi_t> &v) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<N && i <v.size(); i++) block_[i] = v[i];
        }
        /** Initialize the blocks to the values in the argument, b.
            The IndexList, l, specifies the order used to extract
            the indices from b. */
        BlockInfo(const BlockInfo<N> &b, const IndexList &l) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<l.n(); i++) block_[i] = b.block_[l.i(i)];
        }
        template <int NB>
        BlockInfo(const IndexList &l,
                  const BlockInfo<NB> &b,
                  const IndexList &lb) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<l.n(); i++) block_[l.i(i)] = b.block(lb.i(i));
        }
        /** Set all block indices to zero. */
        void zero() { for (int i=0; i<N; i++) block_[i] = 0; }
#ifdef USE_BOUND
        double bound() const { return bound_; }
        void set_bound(double b) const { bound_ = b; }
#endif
        bi_t &block(int i) { return block_[i]; }
        const bi_t &block(int i) const { return block_[i]; }
        /// Compute the size of this block.
        unsigned int size(const Range *indices) const {
          unsigned int r = 1;
          for (int i=0; i<N; i++) r *= indices[i].block_size(block_[i]);
          return r;
        }
        /** Compute the size of a block formed from this block by
            using some subset of the indices given by indexlist. */
        unsigned int subset_size(const Range *indices,
                                 const IndexList &indexlist) const {
          unsigned int r = 1;
          for (int i=0; i<indexlist.n(); i++) {
              int index = indexlist.i(i);
              r *= indices[index].block_size(block_[index]);
            }
          return r;
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il, and another IndexList that gives the index
            mapping into the other BlockInfo, il2. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2, const IndexList &il2) {
          for (int i=0; i<il.n(); i++) {
              block_[il.i(i)] = bi2.block(il2.i(i));
            }
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              block_[il.i(i)] = bi2.block(i);
            }
        }
        /** Return true if blocks are the same as in another BlockInfo,
            bi2, given an IndexList that specifies the index mapping into
            this BlockInfo, il. */
        template <int N2>
        bool equiv_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              if (block_[il.i(i)] != bi2.block(i)) return false;
            }
          return true;
        }
        void print(std::ostream &o=sc::ExEnv::outn()) const;
        void print_block_sizes(const Range *indices,
                               std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<N; i++) {
              if (i) o << " ";
              o << indices[i].block_size(block_[i]);
            }
          o << "}";
        }
        void write(sc::StateOut& so) const {
#ifdef USE_BOUND
          so.put(bound_);
#endif
          for (int i=0; i<N; i++) {
              so.put(int(block_[i]));
            }
        }
        void read(sc::StateIn& si) {
#ifdef USE_BOUND
          si.get(bound_);
#endif
          for (int i=0; i<N; i++) {
              int b;
              si.get(b);
              block_[i] = b;
            }
        }
    };
    template <int N>
    inline std::ostream& operator << (std::ostream&o, const BlockInfo<N> &b)
    {
      b.print(o);
      return o;
    }

    /// \brief Functor for comparing a block's indices.
    template <int N>
    class IndicesLess {
      public:
        bool operator() (const BlockInfo<N>&b1, const BlockInfo<N>&b2) const {
          int b1b = b1.block(0), b2b = b2.block(0);
          if (b1b < b2b) return true;
          if (b1b > b2b) return false;
          for (int i=1; i<N; i++) {
              b1b = b1.block(i); b2b = b2.block(i);
              if (b1b < b2b) return true;
              if (b1b > b2b) return false;
            }
          return false;
        }
        int compare(const BlockInfo<N>&b1, const BlockInfo<N>&b2) const {
          int b1b = b1.block(0), b2b = b2.block(0);
          if (b1b < b2b) return -1;
          if (b1b > b2b) return 1;
          for (int i=1; i<N; i++) {
              b1b = b1.block(i); b2b = b2.block(i);
              if (b1b < b2b) return -1;
              if (b1b > b2b) return 1;
            }
          return 0;
        }
    };

    // Specializations for 0 indices
    template <>
    class BlockInfo<0> {
      private:
#ifdef USE_BOUND
        mutable double bound_;
#endif
      public:
        friend class IndicesLess<0>;
      public:
        BlockInfo() {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        BlockInfo(const std::vector<bi_t> &v) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        BlockInfo(const BlockInfo<2> &b, const IndexList &l) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        template <int NB>
        BlockInfo(const IndexList &l,
                  const BlockInfo<NB> &b,
                  const IndexList &lb) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
#ifdef USE_BOUND
        double bound() const { return bound_; }
        void set_bound(double b) const { bound_ = b; }
#endif
        bi_t &block(int i) {
          throw sc::ProgrammingError("BlockInfo<0>:block: can not be called",
                                     __FILE__, __LINE__);
        }
        const bi_t &block(int i) const {
          throw sc::ProgrammingError("BlockInfo<0>:block: can not be called",
                                     __FILE__, __LINE__);
        }
        /// Compute the size of this block.
        unsigned int size(const Range *indices) const {
          return 1;
        }
        /** Compute the size of a block formed from this block by
            using some subset of the indices given by indexlist. */
        unsigned int subset_size(const Range *indices,
                                 const IndexList &indexlist) const {
          return 1;
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il, and another IndexList that gives the index
            mapping into the other BlockInfo, il2. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2, const IndexList &il2) {
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2) {
        }
        /** Return true if blocks are the same as in another BlockInfo,
            bi2, given an IndexList that specifies the index mapping into
            this BlockInfo, il. */
        template <int N2>
        bool equiv_blocks(const IndexList &il,
                          const BlockInfo<N2> &bi2) {
          return true;
        }
        /** Set all block indices to zero. */
        void zero() {}
        void print(std::ostream &o=sc::ExEnv::outn()) const {
          o << "{}";
        }
        void print_block_sizes(const Range *indices,
                               std::ostream &o=sc::ExEnv::outn()) const {
          o << "{}";
        }
        void write(sc::StateOut& so) const {
#ifdef USE_BOUND
          so.put(bound_);
#endif
        }
        void read(sc::StateIn& si) {
#ifdef USE_BOUND
          si.get(bound_);
#endif
        }
    };
    /// \brief Functor for comparing a block's indices.
    template <>
    class IndicesLess<0> {
      public:
        bool operator() (const BlockInfo<0>&b1, const BlockInfo<0>&b2) const {
          return false;
        }
        int compare(const BlockInfo<0>&b1, const BlockInfo<0>&b2) const {
          return 0;
        }
    };

#if TWO_INDEX_SPECIALIZATIONS
    // Specializations for 2 indices
    template <>
    class BlockInfo<2> {
      private:
#ifdef USE_BOUND
        mutable double bound_;
#endif
      public:
        union {
            sc::sc_uint32_t i;
            sc::sc_uint16_t s[2];
            sc::sc_uint8_t c[4];
        } b_;
        friend class IndicesLess<2>;
      public:
        BlockInfo() {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        BlockInfo(const std::vector<bi_t> &v) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<2 && i<v.size(); i++) b_.s[r2(i)] = v[i];
        }
        BlockInfo(const BlockInfo<2> &b, const IndexList &l) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<l.n(); i++) b_.s[r2(i)] = b.b_.s[r2(l.i(i))];
        }
        template <int NB>
        BlockInfo(const IndexList &l,
                  const BlockInfo<NB> &b,
                  const IndexList &lb) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<l.n(); i++) b_.s[r2(l.i(i))] = b.block(lb.i(i));
        }
#ifdef USE_BOUND
        double bound() const { return bound_; }
        void set_bound(double b) const { bound_ = b; }
#endif
        bi_t &block(int i) { return b_.s[r2(i)]; }
        const bi_t &block(int i) const { return b_.s[r2(i)]; }
        /// Compute the size of this block.
        unsigned int size(const Range *indices) const {
          unsigned int r = 1;
          for (int i=0; i<2; i++) r *= indices[i].block_size(b_.s[r2(i)]);
          return r;
        }
        /** Compute the size of a block formed from this block by
            using some subset of the indices given by indexlist. */
        unsigned int subset_size(const Range *indices,
                                 const IndexList &indexlist) const {
          unsigned int r = 1;
          for (int i=0; i<indexlist.n(); i++) {
              int index = indexlist.i(i);
              r *= indices[index].block_size(b_.s[r2(index)]);
            }
          return r;
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il, and another IndexList that gives the index
            mapping into the other BlockInfo, il2. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2, const IndexList &il2) {
          for (int i=0; i<il.n(); i++) {
              b_.s[r2(il.i(i))] = bi2.block(il2.i(i));
            }
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              b_.s[r2(il.i(i))] = bi2.block(i);
            }
        }
        /** Return true if blocks are the same as in another BlockInfo,
            bi2, given an IndexList that specifies the index mapping into
            this BlockInfo, il. */
        template <int N2>
        bool equiv_blocks(const IndexList &il,
                          const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              if (b_.s[r2(il.i(i))] != bi2.block(i)) return false;
            }
          return true;
        }
        /** Set all block indices to zero. */
        void zero() { b_.i = 0; }
        void print(std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<2; i++) {
              if (i!=0) o << " ";
              o << block(i);
            }
          o << "}";
        }
        void print_block_sizes(const Range *indices,
                               std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<2; i++) {
              if (i) o << " ";
              o << indices[i].block_size(b_.s[r2(i)]);
            }
          o << "}";
        }
        void write(sc::StateOut& so) const {
#ifdef USE_BOUND
          so.put(bound_);
#endif
          for (int i=0; i<2; i++) {
              so.put(int(block(i)));
            }
        }
        void read(sc::StateIn& si) {
#ifdef USE_BOUND
          si.get(bound_);
#endif
          for (int i=0; i<2; i++) {
              int b;
              si.get(b);
              block(i) = b;
            }
        }
    };

    template<> inline
    bool IndicesLess<2>::operator() (const BlockInfo<2>&b1,
                                     const BlockInfo<2>&b2) const
    {
      if (b1.b_.i < b2.b_.i) return true;
      return false;
    }

    template<> inline
    int IndicesLess<2>::compare(const BlockInfo<2>&b1,
                                const BlockInfo<2>&b2) const
    {
      if (b1.b_.i < b2.b_.i) return -1;
      else if (b1.b_.i > b2.b_.i) return 1;
      return 0;
    }
#endif // TWO_INDEX_SPECIALIZATIONS

#if THREE_INDEX_SPECIALIZATIONS
    // Specializations for 3 indices
    template <>
    class BlockInfo<3> {
      private:
#ifdef USE_BOUND
        mutable double bound_;
#endif
      public:
        union {
            sc::sc_uint64_t l;
            sc::sc_uint32_t i[2];
            sc::sc_uint16_t s[4];
            sc::sc_uint8_t c[8];
        } b_;
        friend class IndicesLess<3>;
      public:
        BlockInfo() {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          b_.l = 0;
        }
        BlockInfo(const std::vector<bi_t> &v) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          b_.l = 0;
          for (int i=0; i<3 && i<v.size(); i++) b_.s[r4(i)] = v[i];
        }
        BlockInfo(const BlockInfo<3> &b, const IndexList &l) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          b_.l = 0;
          for (int i=0; i<l.n(); i++) b_.s[r4(i)] = b.b_.s[r4(l.i(i))];
        }
        template <int NB>
        BlockInfo(const IndexList &l,
                  const BlockInfo<NB> &b,
                  const IndexList &lb) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          b_.l = 0;
          for (int i=0; i<l.n(); i++) b_.s[r4(l.i(i))] = b.block(lb.i(i));
        }
#ifdef USE_BOUND
        double bound() const { return bound_; }
        void set_bound(double b) const { bound_ = b; }
#endif
        bi_t &block(int i) { return b_.s[r4(i)]; }
        const bi_t &block(int i) const { return b_.s[r4(i)]; }
        /// Compute the size of this block.
        unsigned int size(const Range *indices) const {
          unsigned int r = 1;
          for (int i=0; i<3; i++) r *= indices[i].block_size(b_.s[r4(i)]);
          return r;
        }
        /** Compute the size of a block formed from this block by
            using some subset of the indices given by indexlist. */
        unsigned int subset_size(const Range *indices,
                                 const IndexList &indexlist) const {
          unsigned int r = 1;
          for (int i=0; i<indexlist.n(); i++) {
              int index = indexlist.i(i);
              r *= indices[index].block_size(b_.s[r4(index)]);
            }
          return r;
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il, and another IndexList that gives the index
            mapping into the other BlockInfo, il2. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2, const IndexList &il2) {
          for (int i=0; i<il.n(); i++) {
              b_.s[r4(il.i(i))] = bi2.block(il2.i(i));
            }
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              b_.s[r4(il.i(i))] = bi2.block(i);
            }
        }
        /** Return true if blocks are the same as in another BlockInfo,
            bi2, given an IndexList that specifies the index mapping into
            this BlockInfo, il. */
        template <int N2>
        bool equiv_blocks(const IndexList &il,
                          const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              if (b_.s[r4(il.i(i))] != bi2.block(i)) return false;
            }
          return true;
        }
        /** Set all block indices to zero. */
        void zero() { b_.l = 0; }
        void print(std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<3; i++) {
              if (i!=0) o << " ";
              o << block(i);
            }
          o << "}";
        }
        void print_block_sizes(const Range *indices,
                               std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<3; i++) {
              if (i) o << " ";
              o << indices[i].block_size(b_.s[r4(i)]);
            }
          o << "}";
        }
        void write(sc::StateOut& so) const {
#ifdef USE_BOUND
          so.put(bound_);
#endif
          for (int i=0; i<3; i++) {
              so.put(int(block(i)));
            }
        }
        void read(sc::StateIn& si) {
#ifdef USE_BOUND
          si.get(bound_);
#endif
          for (int i=0; i<3; i++) {
              int b;
              si.get(b);
              block(i) = b;
            }
        }
    };

    template<> inline
    bool IndicesLess<3>::operator() (const BlockInfo<3>&b1,
                                     const BlockInfo<3>&b2) const
    {
      if (b1.b_.l < b2.b_.l) return true;
      return false;
    }

    template<> inline
    int IndicesLess<3>::compare(const BlockInfo<3>&b1,
                                const BlockInfo<3>&b2) const
    {
      if (b1.b_.l < b2.b_.l) return -1;
      else if (b1.b_.l > b2.b_.l) return 1;
      return 0;
    }
#endif // THREE_INDEX_SPECIALIZATIONS

#if FOUR_INDEX_SPECIALIZATIONS
    // Specializations for 4 indices
    template <>
    class BlockInfo<4> {
      private:
#ifdef USE_BOUND
        mutable double bound_;
#endif
      public:
        union {
            sc::sc_uint64_t l;
            sc::sc_uint32_t i[2];
            sc::sc_uint16_t s[4];
            sc::sc_uint8_t c[8];
        } b_;
        friend class IndicesLess<4>;
      public:
        BlockInfo() {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        BlockInfo(const std::vector<bi_t> &v) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<4 && i<v.size(); i++) b_.s[r4(i)] = v[i];
        }
        BlockInfo(bi_t b0,bi_t b1,bi_t b2,bi_t b3) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          b_.s[0] = b0;
          b_.s[1] = b1;
          b_.s[2] = b2;
          b_.s[3] = b3;
        }
        BlockInfo(const BlockInfo<4> &b, const IndexList &l) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<l.n(); i++) b_.s[r4(i)] = b.b_.s[r4(l.i(i))];
        }
        template <int NB>
        BlockInfo(const IndexList &l,
                  const BlockInfo<NB> &b,
                  const IndexList &lb) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          for (int i=0; i<l.n(); i++) b_.s[r4(l.i(i))] = b.block(lb.i(i));
        }
#ifdef USE_BOUND
        double bound() const { return bound_; }
        void set_bound(double b) const { bound_ = b; }
#endif
        bi_t &block(int i) { return b_.s[r4(i)]; }
        const bi_t &block(int i) const { return b_.s[r4(i)]; }
        /// Compute the size of this block.
        unsigned int size(const Range *indices) const {
          unsigned int r = 1;
          for (int i=0; i<4; i++) r *= indices[i].block_size(b_.s[r4(i)]);
          return r;
        }
        /** Compute the size of a block formed from this block by
            using some subset of the indices given by indexlist. */
        unsigned int subset_size(const Range *indices,
                                 const IndexList &indexlist) const {
          unsigned int r = 1;
          for (int i=0; i<indexlist.n(); i++) {
              int index = indexlist.i(i);
              r *= indices[index].block_size(b_.s[r4(index)]);
            }
          return r;
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il, and another IndexList that gives the index
            mapping into the other BlockInfo, il2. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2, const IndexList &il2) {
          for (int i=0; i<il.n(); i++) {
              b_.s[r4(il.i(i))] = bi2.block(il2.i(i));
            }
        }
        /** Assign blocks to those in another BlockInfo, bi2, given an
            IndexList that specifies the index mapping into this
            BlockInfo, il. */
        template <int N2>
        void assign_blocks(const IndexList &il,
                           const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              b_.s[r4(il.i(i))] = bi2.block(i);
            }
        }
        /** Return true if blocks are the same as in another BlockInfo,
            bi2, given an IndexList that specifies the index mapping into
            this BlockInfo, il. */
        template <int N2>
        bool equiv_blocks(const IndexList &il,
                          const BlockInfo<N2> &bi2) {
          for (int i=0; i<il.n(); i++) {
              if (b_.s[r4(il.i(i))] != bi2.block(i)) return false;
            }
          return true;
        }
        /** Set all block indices to zero. */
        void zero() { b_.l = 0; }
        void print(std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<4; i++) {
              if (i!=0) o << " ";
              o << block(i);
            }
          o << "}";
        }
        void print_block_sizes(const Range *indices,
                               std::ostream &o=sc::ExEnv::outn()) const {
          o << "{";
          for (int i=0; i<4; i++) {
              if (i) o << " ";
              o << indices[i].block_size(b_.s[r4(i)]);
            }
          o << "}";
        }
        void write(sc::StateOut& so) const {
#ifdef USE_BOUND
          so.put(bound_);
#endif
          for (int i=0; i<4; i++) {
              so.put(int(block(i)));
            }
        }
        void read(sc::StateIn& si) {
#ifdef USE_BOUND
          si.get(bound_);
#endif
          for (int i=0; i<4; i++) {
              int b;
              si.get(b);
              block(i) = b;
            }
        }
    };

    template<> inline
    bool IndicesLess<4>::operator() (const BlockInfo<4>&b1,
                                     const BlockInfo<4>&b2) const
    {
      if (b1.b_.l < b2.b_.l) return true;
      return false;
    }

    template<> inline
    int IndicesLess<4>::compare(const BlockInfo<4>&b1,
                                const BlockInfo<4>&b2) const
    {
      if (b1.b_.l < b2.b_.l) return -1;
      else if (b1.b_.l > b2.b_.l) return 1;
      return 0;
    }
#endif // FOUR_INDEX_SPECIALIZATIONS

#ifdef USE_HASH
    // Unary functor for generating hashes from BlockInfo's
    template <int N>
    class BlockInfoHash {
      public:
        size_t operator()(const BlockInfo<N> &b) const {
          size_t r = 0;
          for (int i=0; i<N; i++) r ^= b.block(i);
          return r;
        }
    };
#endif

#ifdef USE_HASH
    // Binary functor for checking if two BlockInfo's are equal
    template <int N>
    class BlockInfoEqual {
      public:
        bool operator()(const BlockInfo<N> &b1,const BlockInfo<N> &b2) const {
          for (int i=0; i<N; i++) if (b1.block(i) != b2.block(i)) return false;
          return true;
        }
    };
#endif

    /** \brief Functor for determining if one IndexList is less than
        another.  Note: this uses the BlockInfo bound info.  As a
        result, when using this as a set or map comparison
        optimization, care must be taken that the bound does not get
        changed, even though the bound is mutable.  */
    template <int N>
    class IndexListLess {
      private:
        IndexList il_;
      public:
        /** Use the IndexList il to sort the indices. */
        IndexListLess(const IndexList &il): il_(il) {}
        /** Use the indices in both IndexList objects, il1 and il2, to
            sort the indices. */
        IndexListLess(const IndexList &il1, const IndexList &il2): il_(il1,il2) {}
        bool operator() (const BlockInfo<N>&b1, const BlockInfo<N>&b2) const {
          if (il_.n() == 0) return false;
          int i0 = il_.i(0);
          int b1b = b1.block(i0), b2b = b2.block(i0);
          if (b1b < b2b) return true;
          if (b1b > b2b) return false;
          for (int l=1; l<il_.n(); l++) {
              int i = il_.i(l);
              b1b = b1.block(i); b2b = b2.block(i);
              if (b1b < b2b) return true;
              if (b1b > b2b) return false;
            }
#ifdef USE_BOUND
          // if indices are equal order according to descending bounds
          if (b1.bound() > b2.bound()) return true;
#endif
          return false;
        }
        int compare(const BlockInfo<N>&b1, const BlockInfo<N>&b2) const {
          if (il_.n() == 0) return 0;
          int i0 = il_.i(0);
          int b1b = b1.block(i0), b2b = b2.block(i0);
          if (b1b < b2b) return -1;
          if (b1b > b2b) return 1;
          for (int l=1; l<il_.n(); l++) {
              int i = il_.i(l);
              b1b = b1.block(i); b2b = b2.block(i);
              if (b1b < b2b) return -1;
              if (b1b > b2b) return 1;
            }
#ifdef USE_BOUND
          // if indices are equal order according to descending bounds
          if (b1.bound() > b2.bound()) return -1;
          else if (b1.bound() < b2.bound()) return 1;
#endif
          return 0;
        }
    };

    /** \brief Data holds the values for each block. */
    class Data: public sc::RefCount {
      private:
//         double *data_;
        size_t ndata_;

        int default_chunksize_;
        // maps n_data available to the pointer and n_data used.
        typedef std::multimap<size_t, std::pair<double *, size_t> > memmap_t;
        memmap_t memmap_;
        typedef std::map<double*,memmap_t::iterator> memitermap_t;
        memitermap_t memitermap_;
      public:
        Data();
        ~Data();
        virtual double *allocate(long size);
        virtual void deallocate(double *);
        double *data() const;
        long ndata() const { return ndata_; }
    };

    /** \brief An Index is used in the symbolic notation for contractions. */
    class Index {
        std::string name_;
        int value_;
      public:
        /** Allocate a named index.  Indices with names can be external
            or internal indices. */
        Index(const std::string &name): name_(name), value_(-1) {}
        /** Allocate an index with a specific value.  The value will
            be fixed during array manipulations.  An Index with a value
            and no name may appear in just a single array, but this
            only makes sense, and is only allowed, when all of the
            blocks for that index are size one.  */
        Index(int value): value_(value) {}
        /** Allocate an Index with a fixed value and a name.  This sort of
            Index can be used to specify fixed external indices.  This is
            the only way to used fixed indices when the block size of the
            fixed index might be greater than one.  The value is the fixed
            block number in this case. */
        Index(const std::string &name, int value): name_(name), value_(value) {}
        /** Return the name of the index.  This is only valid if
            has_name() returns true. */
        const std::string name() const { return name_; }
        /** Return the value of the index.  This is only valid if
            has_value is true. */
        int value() const { return value_; }
        /** Returns true if the index has a value, indicating that
            the index is to be fixed. */
        bool has_value() const { return value_ >= 0; }
        /** Returns true if the indices have non-empty names which are the
            same.  The result is undefined if array operations are done
            when two Index names are the same but their values are
            different. */
        bool symbolically_equivalent(const Index &i) const
            { return name_.size() > 0 && name_ == i.name(); }
        /// Returns true if the indices are equivalent.
        bool operator == (const Index &i) const
            { return name_ == i.name() && value_ == i.value_; }
        /// Set the value of the index.
        void set_value(int value) { value_ = value; }
    };

    template <int N> class Array;
    template <int Nl,int Nr> class ContractProd;
    template <int Nl,int Nr> class ContractUnion;
    /** \brief Represents an array and symbolic indices in a contraction.
        A ContractPart represents the array, a multiplicative factor,
        and indices that occur in symbolic contractions and summations. */
    template <int N>
    class ContractPart {
        Array<N> &array_;
        std::vector<Index> indices_;
        double factor_;
        bool clear_after_use_;
        bool skip_bounds_update_;

        template <int N2, class Op>
        void do_binary_op_with_fixed_indices(double f, const ContractPart<N2> &o,
                                             bool initarray, Op &op) const;

        template <class Op>
        void do_binary_op(double f, const ContractPart<N> &o,
                          bool initarray, Op &op) const;

        template <int Nl,int Nr>
        void doprod(double f, const ContractProd<Nl,Nr> &o,
                   bool initarray) const;

        template <int Nl,int Nr>
        void dounion(const ContractProd<Nl,Nr> &o) const;

      public:
        ContractPart(Array<N> &array);
        ContractPart(Array<N> &array,
                     const Index &i1);
        ContractPart(Array<N> &array,
                     const Index &i1,
                     const Index &i2);
        ContractPart(Array<N> &array,
                     const Index &i1,
                     const Index &i2,
                     const Index &i3);
        ContractPart(Array<N> &array,
                     const Index &i1,
                     const Index &i2,
                     const Index &i3,
                     const Index &i4);
        ContractPart(Array<N> &array,
                     const Index &i1,
                     const Index &i2,
                     const Index &i3,
                     const Index &i4,
                     const Index &i5);
        ContractPart(Array<N> &array,
                     const Index &i1,
                     const Index &i2,
                     const Index &i3,
                     const Index &i4,
                     const Index &i5,
                     const Index &i6);
        void apply_factor(double f);
        double factor() const;
        Array<N>& array() const;
        const Index &index(int i) const;
        bool clear_after_use() const { return clear_after_use_; }

        void operator = (const ContractPart &o) const;
        template <int N2>
        void operator = (const ContractPart<N2> &o) const;
        void operator += (const ContractPart<N> &o) const;
        template <int N2>
        void operator += (const ContractPart<N2> &o) const;
        void operator /= (const ContractPart<N> &o) const;
        template <int N2>
        void operator /= (const ContractPart<N2> &o) const;
        void operator -= (const ContractPart<N> &o) const;
        template <int Nl,int Nr>
        void operator = (const ContractProd<Nl,Nr> &o) const;
        template <int Nl,int Nr>
        void operator += (const ContractProd<Nl,Nr> &o) const;
        template <int Nl,int Nr>
        void operator -= (const ContractProd<Nl,Nr> &o) const;

        /** Add blocks to this corresponding the blocks already allocated
            in o. */
        template <int N2>
        void operator |= (const ContractPart<N2> &o) const;

        template <int Nl,int Nr>
        void operator |= (const ContractProd<Nl,Nr> &o) const;

        /** Instruct the contract routine to clear this array
            immediately after use.  This permits some optimizations to be
            performed. */
        ContractPart<N> operator ~ () const;
        /** Causes the bounds to not be computed for the LHS operand
            for certain operations. */
        ContractPart<N> skip_bounds_update() const;
        /** Extract a value from the array.  All of the indices must be fixed. */
        double value();
    };
    template <int N>
    ContractPart<N> operator *(double, const ContractPart<N> &);

    template <int Nl, int Nr>
    class ContractUnion {
      public:
        ContractPart<Nl> l;
        ContractPart<Nr> r;
        ContractUnion(const ContractPart<Nl> &a1, const ContractPart<Nr> &a2):
          l(a1), r(a2) {}
    };
    template <int Nl, int Nr>
    ContractUnion<Nl,Nr> operator |(const ContractPart<Nl>&l,
                                    const ContractPart<Nr>&r)
    {
      return ContractUnion<Nl,Nr>(l,r);
    }

    template <int N> class BlockDistrib;

    // This checks to see that all elements of vec, vec[i], are
    // in the range [start[i], fence[i]).
    // Used with incr (below).
    template <class T>
    bool
    ready(std::vector<T> &vec,
          const std::vector<T> &start, const std::vector<T> &fence)
    {
      if (vec.size() == 0) return false;
      for (int i=0; i<vec.size(); i++) {
          if (vec[i] < start[i]) return false;
          if (vec[i] >= fence[i]) return false;
        }
      return true;
    }


    // This increments the elements of vec in such a way that all possible
    // vec's with vec[i] in the range [start[i], fence[i]).  are obtained.
    // Used with ready (above).
    template <class T>
    void
    incr(std::vector<T> &vec,
         const std::vector<T> &start, const std::vector<T> &fence)
    {
      if (vec.size() == 0) return;
      for (int i=vec.size()-1; i>0; i--) {
          if (vec[i] >= fence[i]-1) vec[i] = start[i];
          else {
              vec[i]++;
              return;
            }
        }
      vec[0]++;
    }

    /** \brief Implements a block sparse tensor.
        Array maps the BlockInfo to the block's data using the
        given Compare type to sort the blocks. */
    template <int N>
    class Array {
        template <int N1>
        friend
        void remap(typename Array<N1>::cached_blockmap_t& target,
                   const Array<N1> &source,
                   const IndexList &fixed, const BlockInfo<N1> &fixedvals);
        template <int N1>
        friend
        void remap(Array<N1>& target, const Array<N1> &source,
                   const IndexList &il);
      public:
        typedef IndicesLess<N> Compare;
#if USE_STL_MULTIMAP
        typedef typename std::multimap<BlockInfo<N>, double*, Compare > blockmap_t;
        typedef typename std::multimap<BlockInfo<N>, double*, IndexListLess<N> > cached_blockmap_t;
#else
        typedef AVLMMap<BlockInfo<N>, double*, Compare > blockmap_t;
        typedef AVLMMap<BlockInfo<N>, double*, IndexListLess<N> > cached_blockmap_t;
#endif
#ifdef USE_HASH
        typedef typename __gnu_cxx::hash_map<BlockInfo<N>, double*, BlockInfoHash<N>, BlockInfoEqual<N> > blockhash_t;
#endif
      private:
        sc::auto_vec<Range> indices_;
        blockmap_t blocks_;
#ifdef USE_HASH
        blockhash_t block_hash_;
#endif
        sc::Ref<Data> data_;
        bool allocated_;
        std::string name_;
        int debug_;
#ifdef USE_BOUND
        double bound_;
#endif
        double tolerance_;

        std::map<IndexList, cached_blockmap_t*> blockmap_cache_;
        bool use_blockmap_cache_;

        void init_indices(int n) {
          // This routine is used to reset the indices_ to
          // avoid a warning about newing a length 0 array
          // which occurs when N is used.
          if (n) indices_.reset(new Range[n]);
          else indices_.reset(0);
        }

        void basic_init(const std::string &name, double tol) {
          clear();
          use_blockmap_cache_ = true;
          init_indices(N);
          name_ = name;
          debug_ = 0;
          tolerance_ = tol;
          // In the special case that this is a scalar (N == 0),
          // add the block and zero it.
          if (N == 0) {
              add_all_unallocated_blocks();
              zero();
            }
        }
      public:
        ~Array() {
          clear();
        }
        Array(const Array<N> &a) {
          init_indices(N);
          debug_ = 0;
          operator = (a);
        }
        /** Assigns this to 'a', but throws out blocks smaller than
            tol. The tolerance of this will be set to the tolerance of 'a',
            no matter what tol is. */
        void assign_tol(const Array<N> &a,
                        double tol) {
#ifdef USE_BOUND
          bound_ = a.bound_;
#endif
          tolerance_ = a.tolerance_;
          init_blocks(a, tol);
          if (a.allocated_) {
              allocate_blocks();
              for (typename blockmap_t::const_iterator i = a.blocks_.begin();
                   i != a.blocks_.end();
                   i++) {
                  typename blockmap_t::iterator ithis = blocks_.find(i->first);
                  if (ithis == blocks_.end()) continue;
                  memcpy(ithis->second, i->second,
                         sizeof(*ithis->second)*block_size(ithis->first));
                }
              allocated_ = true;
            }
        }
        /** Assigns this to 'a', and keeps all blocks. */
        void assign_all(const Array<N> &a) {
          assign_tol(a,0.0);
        }
        /// Assigns this to 'a', but throws out small blocks
        Array<N> &operator=(const Array<N> &a) {
          assign_tol(a, a.tolerance_);
          return *this;
        }
        /** This does not initialize any indices.  The set_index member
            must be called for all indices before the array is used. */
        Array(const std::string &name = "",
              double tol = DBL_EPSILON) { init(name,tol); }
        void init(const std::string &name = "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
        }
        /** Initialize all N indices to index. */
        Array(const Range &index,
              const std::string & name= "",
              double tol = DBL_EPSILON) { init(index,name,tol); }
        void init(const Range &index,
                  const std::string & name= "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
          if (N < 1)
              throw std::invalid_argument("Array init given too many indices");
          for (int i=0; i<N; i++) set_index(i,index);
        }
        /** Initialize range 0 to i0, and the rest to i1. */
        Array(const Range &i0,const Range &i1,
              const std::string & name= "",
              double tol = DBL_EPSILON) { init(i0,i1,name,tol); }
        void init(const Range &i0,const Range &i1,
                  const std::string & name= "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
          if (N < 2)
              throw std::invalid_argument("Array init given too many indices");
          set_index(0,i0);
          for (int i=1; i<N; i++) set_index(i,i1);
        }
        /** Initialize range 0 to i0, 1 to i1, and the rest to i2. */
        Array(const Range &i0,const Range &i1,const Range &i2,
              const std::string & name= "",
              double tol = DBL_EPSILON) { init(i0,i1,i2,name,tol); }
        void init(const Range &i0,const Range &i1,const Range &i2,
                  const std::string & name= "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
          if (N < 3)
              throw std::invalid_argument("Array init given too many indices");
          set_index(0,i0);
          set_index(1,i1);
          for (int i=2; i<N; i++) set_index(i,i2);
        }
        /** Initialize range 0, 1 and 2 to i0, i1, and i2, and the rest to
            i3. */
        Array(const Range &i0,const Range &i1,
              const Range &i2,const Range &i3,
              const std::string & name= "",
              double tol = DBL_EPSILON) { init(i0,i1,i2,i3,name,tol); }
        void init(const Range &i0,const Range &i1,
                  const Range &i2,const Range &i3,
                  const std::string & name= "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
          if (N < 4)
              throw std::invalid_argument("Array init given too many indices");
          set_index(0,i0);
          set_index(1,i1);
          set_index(2,i2);
          for (int i=3; i<N; i++) set_index(i,i3);
        }
        /** Initialize range 0, 1, 2, 3, and 4 to i0, i1, i2, i3, and i4,
            and the rest to i4. */
        Array(const Range &i0,const Range &i1,
              const Range &i2,const Range &i3,
              const Range &i4,
              const std::string & name= "",
              double tol = DBL_EPSILON) { init(i0,i1,i2,i3,i4,name,tol); }
        void init(const Range &i0,const Range &i1,
                  const Range &i2,const Range &i3,
                  const Range &i4,
                  const std::string & name= "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
          if (N < 5)
              throw std::invalid_argument("Array init given too many indices");
          set_index(0,i0);
          set_index(1,i1);
          set_index(2,i2);
          set_index(3,i3);
          for (int i=4; i<N; i++) set_index(i,i4);
        }
        /** Initialize range 0, 1, 2, 3, 4, and 5 to i0, i1, i2, i3, i4, and i5,
            and the rest to i5. */
        Array(const Range &i0,const Range &i1,
              const Range &i2,const Range &i3,
              const Range &i4,const Range &i5,
              const std::string & name= "",
              double tol = DBL_EPSILON) { init(i0,i1,i2,i3,i4,i5,name,tol); }
        void init(const Range &i0,const Range &i1,
                  const Range &i2,const Range &i3,
                  const Range &i4,const Range &i5,
                  const std::string & name= "",
                  double tol = DBL_EPSILON) {
          basic_init(name,tol);
          if (N < 6)
              throw std::invalid_argument("Array init given too many indices");
          set_index(0,i0);
          set_index(1,i1);
          set_index(2,i2);
          set_index(3,i3);
          set_index(4,i4);
          for (int i=5; i<N; i++) set_index(i,i5);
        }
        /** Initialize each range to the ranges in the given range
            array. */
        Array(const Range *ranges,
              double tol = DBL_EPSILON): data_(new Data), debug_(0),
                                         tolerance_(tol) {
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          set_indices(ranges);
        }
#ifdef USE_BOUND
        double bound() const { return bound_; }
        void set_bound(double b) { bound_ = b; }
#endif
        double tolerance() const { return tolerance_; }
        void set_tolerance(double t) { tolerance_ = t; }
        /** Make this array store nothing. */
        void clear() {
          clear_blockmap_cache();
          allocated_ = false;
          blocks_.clear();
#ifdef USE_HASH
          block_hash_.clear();
#endif
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
          data_ = new Data;
        }
        /// Return the number of indices.
        static int nindex() { return N; }
        /// Return the block map.
        const blockmap_t &blockmap() const { return blocks_; }
#ifdef USE_HASH
        /// Return the block hash.
        const blockhash_t &blockhash() const { return block_hash_; }
#endif
        /** Sets the i'th Range to index. */
        void set_index(int i, const Range &idx) { indices_[i] = idx; }
        /** Initialize each range to the ranges in the given range
            array. */
        void set_indices(const Range *r) {
          for (int i=0; i<N; i++) set_index(i,r[i]);
        }
        /** Gets the i'th Range. */
        const Range &index(int i) const { return indices_.get()[i]; }
        /** Gets the Range array. */
        const Range *indices() const { return indices_.get(); }
        /** Return size of the given block. */
        int block_size(const BlockInfo<N> &bi) const {
          int bs = 1;
          for (int i=0; i<N; i++) bs *= indices_.get()[i].block_size(bi.block(i));
          return bs;
        }
        /** Adds block b.  The data is not allocated.  This cannot be
            called after storage for the data has been allocated with
            allocate blocks. */
        typename blockmap_t::value_type &
        add_unallocated_block(const BlockInfo<N>&b) {
          if (allocated_) {
              throw std::runtime_error(
                  "cannot add unalloced block to alloced array");
            }
          // check the indices
          for (int i=0; i<N; i++) {
              if (b.block(i) >= index(i).nblock()) {
                  throw std::runtime_error(
                      "attempted to add a block with invalid block indices"
                      );
                }
            }
#if USE_STL_MULTIMAP
          typename blockmap_t::iterator bi = blocks_.find(b);
          if (bi != blocks_.end()) return *bi;
          return *blocks_.insert(typename blockmap_t::value_type(b,(double*)0));
#else
          return *blocks_.insert_unique(typename blockmap_t::value_type(b,(double*)0));
#endif
        }
        /** Adds block b.  The data is allocated.  This should not be
            used with add_unallocated_blocks or allocate_blocks.  The
            memory allocated is not initialized.  */
        typename blockmap_t::value_type &
        add_allocated_block(const BlockInfo<N>&b) {
          allocated_ = 1;
          // check the indices
          for (int i=0; i<N; i++) {
              if (b.block(i) >= index(i).nblock()) {
                  throw std::runtime_error(
                      "attempted to add a block with invalid block indices"
                      );
                }
            }
          typename blockmap_t::iterator bi = blocks_.find(b);
          if (bi != blocks_.end()) return *bi;

          double *data = data_->allocate(block_size(b));
          return *blocks_.insert(typename blockmap_t::value_type(b,data));
        }
        /** Adds block b.  The data is allocated.  This should not be
            used with add_unallocated_blocks or allocate_blocks.  The
            data allocated is initialized to zero.  */
        typename blockmap_t::value_type &
        add_zeroed_block(const BlockInfo<N>&b) {
          allocated_ = 1;
          // check the indices
          for (int i=0; i<N; i++) {
              if (b.block(i) >= index(i).nblock()) {
                  throw std::runtime_error(
                      "attempted to add a block with invalid block indices"
                      );
                }
            }

          typename blockmap_t::iterator bi = blocks_.find(b);
          if (bi != blocks_.end()) return *bi;

          int size = block_size(b);
          double *data = data_->allocate(size);
          memset(data, 0, sizeof(double)*size);

          return *blocks_.insert(typename blockmap_t::value_type(b,data));
        }
        /** Removes a block. If the data has been allocated, it is deallocated.
            Does nothing if the block does not exist.
         */
        void remove_block(const BlockInfo<N>&b) {
          typename blockmap_t::iterator bi = blocks_.find(b);
          if (bi != blocks_.end()) {
              data_->deallocate(bi->second);
              blocks_.erase(bi);
            }
        }
        void relocate_block(const BlockInfo<N>&b_old,const BlockInfo<N>&b_new) {
          typename blockmap_t::iterator bi = blocks_.find(b_old);
          if (bi != blocks_.end()) {
              double *data = bi->second;
              blocks_.erase(bi);
              blocks_.insert(std::make_pair(b_new,data));
            }
        }
        /** Adds block b.  The data is not allocated.  This cannot be
            called after storage for the data has been allocated with
            allocate blocks. This assumes that the new block does not
            already exist int the multimap.  The hint gives a suggested
            insertion point. */
        typename blockmap_t::iterator
        add_new_unallocated_block(const typename blockmap_t::iterator &hint,
                                  const BlockInfo<N>&b) {
          if (allocated_) {
              throw std::runtime_error(
                  "cannot add unalloced block to alloced array");
            }
          // check the indices
          for (int i=0; i<N; i++) {
              if (b.block(i) >= index(i).nblock()) {
                  throw std::runtime_error(
                      "attempted to add a block with invalid block indices"
                      );
                }
            }
#if USE_STL_MULTIMAP
          return blocks_.insert(blockmap_t::value_type(b,(double*)0));
#else
          return blocks_.insert_new(hint, typename blockmap_t::value_type(b,(double*)0));
#endif
        }
        /** Adds block b.  The data is not allocated.  This cannot be
            called after storage for the data has been allocated with
            allocate blocks. This assumes that the new block does not
            already exist int the multimap.  */
        typename blockmap_t::iterator
        add_new_unallocated_block(const BlockInfo<N>&b) {
          if (allocated_) {
              throw std::runtime_error(
                  "cannot add unalloced block to alloced array");
            }
          // check the indices
          for (int i=0; i<N; i++) {
              if (b.block(i) > index(i).nblock()) {
                  throw std::runtime_error(
                      "attempted to add a block with invalid block indices"
                      );
                }
            }
#if USE_STL_MULTIMAP
          return blocks_.insert(blockmap_t::value_type(b,(double*)0));
#else
          return blocks_.insert_new(typename blockmap_t::value_type(b,(double*)0));
#endif
        }
        /** Adds all possible blocks to give a dense array.  The data is
            not allocated.  This cannot be called after storage for the
            data has been allocated with allocate blocks.  */
        void
        add_all_unallocated_blocks() {

          if (allocated_) {
              throw std::runtime_error(
                  "cannot add unalloced block to alloced array");
            }

          clear();

          if (N == 0) {
              BlockInfo<N> bi;
              add_new_unallocated_block(bi);
              return;
            }

          std::vector<sma2::bi_t> b_indices(N);
          std::vector<sma2::bi_t> b_start(N);
          std::vector<sma2::bi_t> b_fence(N);

          std::fill(b_start.begin(), b_start.end(), 0);
          for (int i=0; i<N; i++) b_fence[i] = index(i).nblock();

          for (b_indices=b_start;
               ready(b_indices, b_start, b_fence);
               incr(b_indices, b_start, b_fence)) {
              sma2::BlockInfo<N> bi(b_indices);
              add_new_unallocated_block(bi);
            }
        }
        /** Returns an initial hint. Basically, it is just used to
           initialize a blockmap::iterator to a valid value. */
        typename blockmap_t::iterator initial_hint() {
          return blocks_.end();
        }
        /** Returns the number of blocks. */
        int n_block() const {
          return blocks_.size();
        }
        /** Returns the number of elements contained in blocks. */
        size_t n_element() const {
          size_t n = 0;
          for (typename blockmap_t::const_iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              n += block_size(i->first);
            }
          return n;
        }
        /** Returns the number of elements contained in blocks. */
        size_t n_element_allocated() const {
          if (data_ == 0) return 0;
          return data_->ndata();
        }
        std::string name() const { return name_; }
        double max_n_block() {
          double n_block_max = 1;
          for (int i=0; i<N; i++) {
              n_block_max *= index(i).nblock();
            }
          return n_block_max;
        }
        double max_n_element() {
          double n_element_max = 1;
          for (int i=0; i<N; i++) {
              n_element_max *= index(i).nindex();
            }
          return n_element_max;
        }
        /** Assigns the given value to all elements in the array. */
        void assign(double val) {
          allocate_blocks();
          for (typename blockmap_t::iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              int n = block_size(i->first);
              double *d = i->second;
              for (int j=0; j<n; j++) d[j] = val;
#ifdef USE_BOUND
              i->first.set_bound(val);
#endif
            }
#ifdef USE_BOUND
          bound_ = val;
#endif
        }
        /** Zeros the array. */
        void zero() {
          assign(0.0);
        }
        /** Computes the bound for each block. */
        void compute_bounds() {
          allocate_blocks();
#ifdef USE_BOUND
          bound_ = 0.0;
          for (typename blockmap_t::iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              int n = block_size(i->first);
              double *d = i->second;
              double maxval = 0.0;
              for (int j=0; j<n; j++) {
                  double val = fabs(d[j]);
                  if (val > maxval) maxval = val;
                }
              i->first.set_bound(maxval);
              if (maxval > bound_) bound_ = maxval;
            }
#endif
        }
        /** Allocate storage for all blocks if storage has not already been
            allocated. */
        void allocate_blocks() {
          if (allocated_) return;
          size_t nblock = 0;
          long n = n_element();
          // clj debug
          //if (n > 1000000) {
          //    sc::Ref<sc::MessageGrp> msg
          //        = sc::MessageGrp::get_default_messagegrp();
          //    int me = msg->me();
          //    sc::ExEnv::outn() << me << ": " << name()
          //              << " allocating " << n << std::endl;
          //  }
          // clj end debug
          double *dat = data_->allocate(n);
          for (typename blockmap_t::iterator i = blocks_.begin();
               i != blocks_.end();
               i++,nblock++) {
              if (i->second == 0) {
                  int n = block_size(i->first);
                  i->second = dat;
                  dat += n;
                }
              else {
                  throw std::runtime_error("allocate_blocks: duplicate alloc");
                }
            }
#ifdef USE_HASH
          block_hash_.resize(nblock);
          for (typename blockmap_t::iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              block_hash_.insert(*i);
            }
#endif
          allocated_ = true;
        }
        /** Deallocate storage for all blocks. */
        void deallocate_blocks() {
          data_ = new Data;
          for (typename blockmap_t::iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              i->second = 0;
            }
          allocated_ = false;
#ifdef USE_BOUND
          bound_ = DBL_MAX;
#endif
        }
        ContractPart<N> operator()() {
          return ContractPart<N>(*this);
        }
        ContractPart<N> operator()(const Index &i1) {
          return ContractPart<N>(*this,i1);
        }
        ContractPart<N> operator()(const Index &i1, const Index &i2) {
          return ContractPart<N>(*this,i1,i2);
        }
        ContractPart<N> operator()(const std::string &i1,
                                   const std::string &i2) {
          return ContractPart<N>(*this,i1,i2);
        }
        ContractPart<N> operator()(const Index &i1, const Index &i2, const Index &i3) {
          return ContractPart<N>(*this,i1,i2,i3);
        }
        ContractPart<N> operator()(const std::string &i1,
                                   const std::string &i2,
                                   const std::string &i3) {
          return ContractPart<N>(*this,i1,i2,i3);
        }
        ContractPart<N> operator()(const Index &i1, const Index &i2,
                                   const Index &i3, const Index &i4) {
          return ContractPart<N>(*this,i1,i2,i3,i4);
        }
        ContractPart<N> operator()(const std::string &i1,
                                   const std::string &i2,
                                   const std::string &i3,
                                   const std::string &i4) {
          return ContractPart<N>(*this,i1,i2,i3,i4);
        }
        ContractPart<N> operator()(const Index &i1, const Index &i2,
                                   const Index &i3, const Index &i4,
                                   const Index &i5) {
          return ContractPart<N>(*this,i1,i2,i3,i4,i5);
        }
        ContractPart<N> operator()(const std::string &i1,
                                   const std::string &i2,
                                   const std::string &i3,
                                   const std::string &i4,
                                   const std::string &i5) {
          return ContractPart<N>(*this,i1,i2,i3,i4,i5);
        }
        ContractPart<N> operator()(const Index &i1, const Index &i2,
                                   const Index &i3, const Index &i4,
                                   const Index &i5, const Index &i6) {
          return ContractPart<N>(*this,i1,i2,i3,i4,i5,i6);
        }
        ContractPart<N> operator()(const std::string &i1,
                                   const std::string &i2,
                                   const std::string &i3,
                                   const std::string &i4,
                                   const std::string &i5,
                                   const std::string &i6) {
          return ContractPart<N>(*this,i1,i2,i3,i4,i5,i6);
        }
        /** Multiplication by a scalar */
        void operator *= (double f) {
          if (f == 0.0) { clear(); compute_bounds(); return; }
          for (typename blockmap_t::const_iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              double fabs_f = fabs(f);
              int n = block_size(i->first);
              double *data = i->second;
              for (int j=0; j<n; j++) data[j] *= f;
#ifdef USE_BOUND
              i->first.set_bound(fabs_f * i->first.bound());
#endif
            }
        }
        /** Initialize the stored blocks to be the same as those in a. */
        void init_blocks(const Array<N> &a, double tol) {
          clear();
          set_indices(a.indices());
          for (typename blockmap_t::const_iterator i = a.blocks_.begin();
               i != a.blocks_.end();
               i++) {
              if (a.block_max_abs(i) > tol) {
                  add_unallocated_block(i->first);
                }
            }
        }
        /** Initialize the stored blocks to be the same as those in a. */
        void init_blocks(const Array<N> &a) {
          init_blocks(a, tolerance());
        }
        /** Find maximum absolute value of elements in a block. */
        double block_max_abs(typename blockmap_t::const_iterator &b) const {
          int s = block_size(b->first);
          const double *dat = b->second;
          double maxabs = 0.0;
          for (int i=0; i<s; i++) {
              double v = fabs(dat[i]);
              if (v > maxabs) maxabs = v;
            }
          return maxabs;
        }
        /** Find maximum absolute value of elements */
        double max_abs_element() {
          double max_abs = 0.0;
          for (typename blockmap_t::const_iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
            max_abs = std::max(max_abs, this->block_max_abs(i));
          }
          return max_abs;
        }

        /** Print the array. */
        void print_local(std::ostream&o=sc::ExEnv::outn()) const;
        void print(const sc::Ref<sc::MessageGrp> &grp = 0,
                   bool distributed = false, std::ostream&o=sc::ExEnv::outn()) const;

        /** Set/get the debug level. */
        int &debug() { return debug_; }
        const int &debug() const { return debug_; }

        void read(sc::StateIn&si) {
          clear();
          for (int i=0; i<N; i++) {
              indices_[i].read(si);
              indices_[i].print();
            }
          int allocatedval;
          si.get(allocatedval);
          char *str;
          si.getstring(str);
          name_ = str;
          delete[] str;
#ifdef USE_BOUND
          si.get(bound_);
          sc::ExEnv::outn() << "bound = " << bound_ << std::endl;
#endif
          si.get(tolerance_);
          sc::ExEnv::outn() << "tolerance = " << tolerance_ << std::endl;
          int nblock;
          si.get(nblock);
          BlockInfo<N> bi;
          for (int i=0; i<nblock; i++) {
              bi.read(si);
              add_unallocated_block(bi);
            }
          if (allocatedval) {
              allocate_blocks();
              si.get_array_double(data_->data(), data_->ndata());
            }
        }
        void write(sc::StateOut&so) const {
          for (int i=0; i<N; i++) indices_[i].write(so);
          int bval = allocated_;
          so.put(bval);
          so.putstring(name_.c_str());
#ifdef USE_BOUND
          so.put(bound_);
#endif
          so.put(tolerance_);
          int blocks_size = int(blocks_.size());
          if (blocks_.size() != blocks_size) {
            throw std::runtime_error("Array::write: blocks truncated");
          }
          so.put(blocks_size);
          for (typename blockmap_t::const_iterator i = blocks_.begin();
               i != blocks_.end();
               i++) {
              i->first.write(so);
            }
          if (allocated_) so.put_array_double(data_->data(), data_->ndata());
        }
        /// \brief Performs a reduce-broadcast on the data.
        /// The array must have the same block structure on all processes.
        void parallel_accumulate(const sc::Ref<sc::MessageGrp> &grp);
        /// \brief Makes the block distribution the same on all processes.
        /// If any processes hold a block, then all will after this call.
        /// This does not send any data, only block structure is sent.
        void parallel_union(const sc::Ref<sc::MessageGrp> &grp);
        /// Convert a distributed array to a replicated array.
        /// The result is stored in this.
        void replicated_from_distributed(const sc::Ref<sc::MessageGrp>&,
                                         const Array<N> &);
        /// Redistribute an array.
        /// The result is stored in this.
        void distributed_from_distributed(const sc::Ref<sc::MessageGrp>&,
                                          const BlockDistrib<N> &,
                                          Array<N> &,
                                          bool clear_source_array = false,
                                          bool ignore_block_distrib_throws = false);

        /// If this is a scalar (N==0) return the value of the scalar.
        double value() {
          if (N != 0) throw sc::ProgrammingError("Array:value: N!=0",
                                                 __FILE__, __LINE__);
          return blockmap().begin()->second[0];
        }
        /** Get rid of cached blockmaps. */
        void clear_blockmap_cache() {
          for (typename std::map<IndexList,cached_blockmap_t*>::iterator
                   iter = blockmap_cache_.begin();
               iter != blockmap_cache_.end();
               iter++) {
              delete iter->second;
            }
          blockmap_cache_.clear();
        }
        /// \brief Used to indicate whether or not a blockmap cache is to be used.
        /// Blockmap caches store alternative blockmaps that are used in
        /// the contract routine. If a required blockmap is not available, then
        /// it is computed. It is also added to the blockmap cache if this is
        /// true.
        void set_use_blockmap_cache(bool ubmc) {
          use_blockmap_cache_ = ubmc;
        }
        /// Return true if blockmap caches are in use.
        bool use_blockmap_cache() const {
          return use_blockmap_cache_;
        }
        /// Return true if the needed blockmap cache entry exists.
        bool blockmap_cache_entry_exists(const IndexList &il) {
          return blockmap_cache_.find(il) != blockmap_cache_.end();
        }
        /// Return a cached blockmap.
        cached_blockmap_t &blockmap_cache_entry(const IndexList &il) {
          cached_blockmap_t *&entry = blockmap_cache_[il];
          if (entry == 0) {
              entry = new cached_blockmap_t(il);
              IndexList nullil;
              BlockInfo<N> nullbi;
              remap(*entry,*this,nullil,nullbi);
            }
          return *entry;
        }
    };
    template <int N>
    inline std::ostream& operator << (std::ostream&o, const Array<N> &a)
    {
      a.print(o);
      return o;
    }

    inline int offset2(int i, int ni, int j, int nj) {
      return i*nj+j;
    }

    template <int N> double scalar_contract(Array<N> &c, Array<N> &a, const IndexList &alist);

    /** \brief Represents a pairs of contracted array and their
        symbolic indices.
     */
    template <int Nl, int Nr>
    class ContractProd {
      public:
        sc::Ref<sc::RegionTimer> regtimer;
        ContractPart<Nl> l;
        ContractPart<Nr> r;
        void apply_factor(double f) { l.apply_factor(f); }
        ContractProd(const ContractPart<Nl> &a1, const ContractPart<Nr> &a2):
          l(a1), r(a2) {}
        operator double() {
          std::vector<int> ivec;
          for (int i=0; i<Nr; i++) {
              for (int j=0; j<Nl; j++) {
                  if (r.index(i) == l.index(j)) {
                      ivec.push_back(j);
                      break;
                    }
                }
            }
          IndexList il(ivec);
          return l.factor() * r.factor()
              * sma2::scalar_contract(l.array(), r.array(), il);
        }
        ContractProd timer(const sc::Ref<sc::RegionTimer> &t) {
          ContractProd<Nl,Nr> ret(*this);
          ret.regtimer = t;
          return ret;
        }
    };
    template <int Nl, int Nr>
    ContractProd<Nl,Nr> operator *(const ContractPart<Nl>&l,
                                   const ContractPart<Nr>&r)
    {
      return ContractProd<Nl,Nr>(l,r);
    }
    template <int Nl, int Nr>
    ContractProd<Nl,Nr> operator *(double f, const ContractProd<Nl,Nr> &prod)
    {
      ContractProd<Nl,Nr> result(prod);
      result.apply_factor(f);
      return result;
    }

    /** Remap the data using a different index ordering.  This does not
        repack the data, so the remapped indices must have a block size of
        one.  This creates a reference to the original array's data, so the
        data cannot be changed. */
    template <int N>
    void
    remap(Array<N>& target, const Array<N> &source,
          const IndexList &il)
    {
      for (int i=0; i<N; i++) target.set_index(i,source.index(i));
      const typename Array<N>::blockmap_t &sourceblocks
          = source.blockmap();
      target.blocks_.clear();
      for (typename Array<N>::blockmap_t::const_iterator i = sourceblocks.begin();
           i != sourceblocks.end();
           i++) {
          const BlockInfo<N> &orig_bi(i->first);
          BlockInfo<N> new_bi(orig_bi,il);
          target.blocks_.insert(std::pair<BlockInfo<N>,double*>(new_bi,
                                                                     i->second));
        }
      target.data_ = source.data_;
    }

    /** Remap the data using a different comparision operation.  Only
        blocks with the given fixed index are put into the new map.  The
        fixed indices must be the first indices in the source array.  This
        creates a reference to the original array's data, so the data
        cannot be changed. */
    template <int N>
    void
    remap(typename Array<N>::cached_blockmap_t& target,
          const Array<N> &source,
          const IndexList &fixed, const BlockInfo<N> &fixedvals)
    {
      if (!fixed.is_identity_permutation())
          throw std::invalid_argument("remap requires fixed indices first");
      const typename Array<N>::blockmap_t &sourceblocks
          = source.blockmap();
      target.clear();
      typename Array<N>::blockmap_t::const_iterator begin, end;
      if (fixed.n() == 0) {
          begin = sourceblocks.begin();
          end = sourceblocks.end();
        }
      else {
          BlockInfo<N> sbi;

          sbi.zero();
          sbi.assign_blocks(fixed, fixedvals);
#ifdef USE_BOUND
          sbi.set_bound(DBL_MAX);
#endif
          begin = sourceblocks.lower_bound(sbi);

          for (int i=0; i<N; i++) sbi.block(i) = source.index(i).nblock();
          sbi.assign_blocks(fixed, fixedvals);
#ifdef USE_BOUND
          sbi.set_bound(0.0);
#endif
          end = sourceblocks.upper_bound(sbi);
        }
      for (typename Array<N>::blockmap_t::const_iterator i = begin;
           i != end;
           i++) {
          target.insert(*i);
        }
    }

    /** Puts the diagonal elements of the array in diag. */
    void extract_diagonal(Array<2> &array, std::vector<double> &diag);

    /** Scales the diagonal elements of array. */
    void scale_diagonal(Array<2> &array, double factor);

    /** Scales the elements of the array with an MP2 denominator. V must have
        an operator() member to access to the elements. */
    template <int N, class V>
    void
    apply_denominator(Array<N> &array, double denominator_offset,
                      const std::vector<V> &eigenvalues);

    /** \brief BlockIter loops through the all the indices within a
        block. */
    template <int N>
    class BlockIter {
        int bs_[N];
        int i_[N];
        bool ready_;
        void init(const Range *index,
                  const BlockInfo<N> &bi,
                  const IndexList &il) {
          for (int i=0; i<N; i++) {
              bs_[i] = index[il.i(i)].block_size(bi.block(il.i(i)));
            }
        }
      public:
        BlockIter(const Range *index,
                  const BlockInfo<N> &bi) {
          init(index,bi,IndexList::identity(N));
        }
        BlockIter(const Range *index,
                  const BlockInfo<N> &bi,
                  const IndexList &il) {
          init(index,bi,il);
        }
        /// Initialize the iterator.
        void start() { ready_ = true; for (int i=0; i<N; i++) i_[i] = 0; }
        /// Returns true if the iterator is valid.
        bool ready() const { return ready_; }
        /// Increments the iterator.
        void operator++(int) {
          for (int i=N-1; i>=0; i--) {
              i_[i]++;
              if (i_[i] == bs_[i]) {
                  i_[i] = 0;
                }
              else return;
            }
          ready_ = false;
        }
        int offset() {
          int r = 0;
          for (int i=0; i<N; i++) {
              r = r * bs_[i] + i_[i];
            }
          return r;
        }
        int subset_offset(const IndexList &subset_indexlist) {
          int r = 0;
          for (int i=0; i<subset_indexlist.n(); i++) {
              r = r * bs_[subset_indexlist.i(i)]
                + i_[subset_indexlist.i(i)];
            }
          return r;
        }
        /// Returns the size of the given block.
        int block_size(int i) const { return bs_[i]; }
        int &index(int i) { return i_[i]; }
        const int &index(int i) const { return i_[i]; }
    };

    /** \brief Blocksize == 0 specialization of BlockIter. */
    template <>
    class BlockIter<0> {
        bool ready_;
      public:
        BlockIter(const Range *index,
                  const BlockInfo<0> &bi) {}
        BlockIter(const Range *index,
                  const BlockInfo<0> &bi,
                  const IndexList &il) {}
        void start() { ready_ = true; }
        bool ready() const { return ready_; }
        void operator++(int) {
          ready_ = false;
        }
        int offset() {
          return 0;
        }
        int subset_offset(const IndexList &subset_indexlist) {
          return 0;
        }
        int block_size(int i) const {
          throw sc::ProgrammingError("BlockIter<0>:block_size: can not be called");
        }
        int &index(int i) {
          throw sc::ProgrammingError("BlockIter<0>:index: can not be called");
        }
        const int &index(int i) const {
          throw sc::ProgrammingError("BlockIter<0>:index: can not be called");
        }
    };

    /** The sum function computes the sum of two arrays. */
    template <int N, class Op>
    void binary_op(Array<N> &c,
                   double alpha, const Array<N> &a, const IndexList &alist,
                   bool skip_bounds_update, Op &op)
    {
      // Consistency checks
      if (alist.n() != N) {
          throw std::invalid_argument("sma2::binary_op: # of indices inconsistent");
        }
//        sc::ExEnv::outn() << "a's indices on c" << std::endl;
//        for (int i=0; i<N; i++) {
//            sc::ExEnv::outn() << " " << alist.i(i);
//          }
//        sc::ExEnv::outn() << std::endl;
      for (int i=0; i<N; i++) {
          if (c.index(alist.i(i)) != a.index(i))
              throw std::invalid_argument("sma2::binary_op: indices don't agree");
        }

      bool same_index_order = alist.is_identity();

      // if c has fewer blocks then it drives the summation
      // otherwise a drives
      if (c.n_block() < a.n_block()) {
          const typename Array<N>::blockmap_t &amap = a.blockmap();
          const typename Array<N>::blockmap_t &cmap = c.blockmap();
          IndexList clist = alist.reverse_mapping();
          for (typename Array<N>::blockmap_t::const_iterator citer = cmap.begin();
               citer != cmap.end();
               citer++) {
              const BlockInfo<N> &cbi = citer->first;
              BlockInfo<N> abi(cbi,alist);
              //sc::ExEnv::outn() << "summing into " << cbi << std::endl;
              typename Array<N>::blockmap_t::const_iterator ablock
                  = amap.find(abi);
              if (ablock == amap.end()) continue;
              double *cdata = citer->second;
              double *adata = ablock->second;
              if (same_index_order) {
                  int sz = c.block_size(cbi);
                  for (int i=0; i<sz; i++) op(cdata[i],alpha * adata[i]);
                }
              else {
                  BlockIter<N> cbiter(c.indices(),cbi);
                  int coff = 0;
                  for (cbiter.start(); cbiter.ready(); cbiter++,coff++) {
                      op(cdata[coff],
                         alpha * adata[cbiter.subset_offset(alist)]);
                    }
                }
            }
        }
      else {
          const typename Array<N>::blockmap_t &amap = a.blockmap();
          const typename Array<N>::blockmap_t &cmap = c.blockmap();
          IndexList clist = alist.reverse_mapping();
          for (typename Array<N>::blockmap_t::const_iterator aiter = amap.begin();
               aiter != amap.end();
               aiter++) {
              const BlockInfo<N> &abi = aiter->first;
              BlockInfo<N> cbi(abi,clist);
              typename Array<N>::blockmap_t::const_iterator cblock
                  = cmap.find(cbi);
              if (cblock == cmap.end()) continue;
              double *adata = aiter->second;
              double *cdata = cblock->second;
              if (same_index_order) {
                  int sz = a.block_size(abi);
                  for (int i=0; i<sz; i++) op(cdata[i], alpha * adata[i]);
                }
              else {
                  BlockIter<N> cbiter(c.indices(),cbi);
                  int coff = 0;
                  for (cbiter.start(); cbiter.ready(); cbiter++,coff++) {
                      op(cdata[coff],
                         alpha * adata[cbiter.subset_offset(alist)]);
                    }
                }
            }
        }
      if (!skip_bounds_update) c.compute_bounds();
    }

    /** Packs a matrix into a sma2::Array.  Class T must have a double
        operator()(int,int) member.  The nonzero blocks of the array must
        have already been specified. */
    template <class T>
    void pack_matrix_into_array(const T &m, Array<2> &smaa) {
      smaa.allocate_blocks();
      const Array<2>::blockmap_t &Abm = smaa.blockmap();
      for (Array<2>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<2> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          int off1 = smaa.index(1).block_offset(Abi.block(1));
          BlockIter<2> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              dat[j.offset()] = m(off0+j.index(0), off1+j.index(1));
            }
        }
      smaa.compute_bounds();
    }

    /** Packs a vector into a sma2::Array.  Class T must have a double
        operator()(int) member.  The nonzero blocks of the array must have
        already been specified. */
    template <class T>
    void pack_vector_into_array(const T &m, Array<1> &smaa) {
      smaa.allocate_blocks();
      const Array<1>::blockmap_t &Abm = smaa.blockmap();
      for (Array<1>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<1> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          BlockIter<1> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              dat[j.offset()] = m(off0+j.index(0));
            }
        }
      smaa.compute_bounds();
    }

    /** Packs a matrix into a sma2::Array.  Class T must have a double
        operator()(int,int) member.  The nonzero blocks of the array must
        have already been specified. The matrix does not have the full
        dimensions represented by the ranges of the array.  The index_map
        specifies how to convert the matrix's indices to the array's
        indices.  No checking is done, the Array's blocks must refer to
        only valid elements of the matrix.  */
    template <class T>
    void pack_submatrix_into_array(const T &m, Array<2> &smaa,
                                   const std::vector<int> &index_map_i,
                                   const std::vector<int> &index_map_j) {

      std::vector<int> reverse_map_i(smaa.index(0).nindex());
      std::vector<int> reverse_map_j(smaa.index(1).nindex());
      std::fill(reverse_map_i.begin(), reverse_map_i.end(), -1);
      std::fill(reverse_map_j.begin(), reverse_map_j.end(), -1);
      for (int i=0; i<index_map_i.size(); i++) {
          reverse_map_i[index_map_i[i]] = i;
        }
      for (int i=0; i<index_map_j.size(); i++) {
          reverse_map_j[index_map_j[i]] = i;
        }

      smaa.allocate_blocks();
      const Array<2>::blockmap_t &Abm = smaa.blockmap();
      for (Array<2>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<2> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          int off1 = smaa.index(1).block_offset(Abi.block(1));
          BlockIter<2> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              int m_i = reverse_map_i[off0+j.index(0)];
              int m_j = reverse_map_j[off1+j.index(1)];
              if (m_i == -1 || m_j == -1) continue;
              dat[j.offset()] = m(m_i, m_j);
            }
        }
      smaa.compute_bounds();
    }

    /** Packs a vector into a sma2::Array.  Class T must have a double
        operator()(int) member.  The nonzero blocks of the array must
        have been specified in advance. The vector does not have the full
        dimensions represented by the ranges of the array.  The index_map
        specifies how to convert the vector's indices to the array's
        indices.  No checking is done, the Array's blocks must refer to
        only valid elements of the matrix.  */
    template <class T>
    void pack_subvector_into_array(const T &m, Array<1> &smaa,
                                   const std::vector<int> &index_map) {

      std::vector<int> reverse_map(smaa.index(0).nindex());
      std::fill(reverse_map.begin(), reverse_map.end(), -1);
      for (int i=0; i<index_map.size(); i++) {
          reverse_map[index_map[i]] = i;
        }

      smaa.allocate_blocks();
      const Array<1>::blockmap_t &Abm = smaa.blockmap();
      for (Array<1>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<1> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          BlockIter<1> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              int m_i = reverse_map[off0+j.index(0)];
              if (m_i == -1) continue;
              dat[j.offset()] = m(m_i);
            }
        }
      smaa.compute_bounds();
    }

    /** Packs a dense matrix into a sma2::Array.  Class T must have a
        double operator()(int,int) member.  The matrix m corresponds to a
        subset of a full matrix.  The blocks that make up this subset are
        given in the row_blcoks and col_blocks arguments. */
    template <class T>
    void pack_dense_matrix_into_array(const T &m, Array<2> &smaa,
                                      const std::set<int> &row_blocks,
                                      const std::set<int> &col_blocks) {
      // allocate the nonzero blocks of smaa
      smaa.clear();
      for (std::set<int>::const_iterator it1 = row_blocks.begin();
           it1 != row_blocks.end(); it1++) {
          for (std::set<int>::const_iterator it2 = col_blocks.begin();
               it2 != col_blocks.end(); it2++) {
              BlockInfo<2> bi;
              bi.block(0) = *it1; bi.block(1) = *it2;
              smaa.add_unallocated_block(bi);
            }
        }
      smaa.allocate_blocks();

      // initialize the offset arrays
      int off=0;
      std::vector<int> row_offset(smaa.index(0).nblock());
      for (std::set<int>::const_iterator i = row_blocks.begin();
           i!=row_blocks.end(); i++) {
          row_offset[*i] = off;
          off+=smaa.index(0).block_size(*i);
        }
      off=0;
      std::vector<int> col_offset(smaa.index(1).nblock());
      for (std::set<int>::const_iterator i = col_blocks.begin();
           i!=col_blocks.end(); i++) {
          col_offset[*i] = off;
          off+=smaa.index(1).block_size(*i);
        }

      // place the data in the array
      const Array<2>::blockmap_t &Abm = smaa.blockmap();
      for (Array<2>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<2> &Abi = i->first;
          int off0 = row_offset[Abi.block(0)];
          int off1 = col_offset[Abi.block(1)];
          BlockIter<2> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              dat[j.offset()] = m(off0+j.index(0), off1+j.index(1));
            }
        }
      smaa.compute_bounds();
    }

    /** Packs a matrix into an empty sma2::Array.  Class T must have a
        double operator()(int,int) member.  The nonzero blocks of the array
        will be determined. */
    template <class T>
    void pack_matrix_into_empty_array(const T &m, Array<2> &smaa,
                                double bound) {
      smaa.clear();
      const Range &index0 = smaa.index(0);
      const Range &index1 = smaa.index(1);
      for (int i0=0; i0 < index0.nblock(); i0++) {
          int bs0 = index0.block_size(i0);
          int bo0 = index0.block_offset(i0);
          for (int i1=0; i1 < index1.nblock(); i1++) {
              int bs1 = index1.block_size(i1);
              int bo1 = index1.block_offset(i1);
              double maxval = 0.0;
              for (int j0=0; j0<bs0; j0++) {
                  for (int j1=0; j1<bs1; j1++) {
                      double val = fabs(m(bo0+j0,bo1+j1));
                      if (val > maxval) maxval = val;
                    }
                }
              if (maxval >= bound) {
                  BlockInfo<2> bi;
                  bi.block(0) = i0;
                  bi.block(1) = i1;
                  double *data = smaa.add_unallocated_block(bi).second;
                }
            }
        }
      pack_matrix_into_array(m,smaa);
    }

    /** Packs a vector into an empty sma2::Array.  Class T must have a
        double operator()(int) member and must support dim().n().  The
        nonzero blocks of the array will be determined. */
    template <class T>
    void pack_vector_into_empty_array(const T &m, Array<1> &smaa,
                                      double bound) {
      if (m.dim().n() != smaa.index(0).nindex()) {
          throw std::runtime_error("pack_vector_into_empty_array: dims don't match");
        }
      smaa.clear();
      const Range &index0 = smaa.index(0);
      for (int i0=0; i0 < index0.nblock(); i0++) {
          int bs0 = index0.block_size(i0);
          int bo0 = index0.block_offset(i0);
          double maxval = 0.0;
          for (int j0=0; j0<bs0; j0++) {
              double val = fabs(m(bo0+j0));
              if (val > maxval) maxval = val;
            }
          if (maxval >= bound) {
              BlockInfo<1> bi;
              bi.block(0) = i0;
              double *data = smaa.add_unallocated_block(bi).second;
            }
        }
      pack_vector_into_array(m,smaa);
    }

    /** Packs a matrix into an empty sma2::Array.  Class T must have a
        double operator()(int,int) member.  The nonzero blocks of the array
        will be determined. The matrix does not have the full dimensions
        represented by the ranges of the array.  The index_map specifies
        how to convert the matrix's indices to the array's indices. */
    template <class T>
    void pack_submatrix_into_empty_array(const T &m, Array<2> &smaa,
                                         double bound,
                                         const std::vector<int> &index_map_i,
                                         const std::vector<int> &index_map_j) {
      smaa.clear();
      const Range &index0 = smaa.index(0);
      const Range &index1 = smaa.index(1);
      std::vector<int> reverse_map_i(index0.nindex());
      std::vector<int> reverse_map_j(index1.nindex());
      std::fill(reverse_map_i.begin(), reverse_map_i.end(), -1);
      std::fill(reverse_map_j.begin(), reverse_map_j.end(), -1);
      for (int i=0; i<index_map_i.size(); i++) {
          reverse_map_i[index_map_i[i]] = i;
        }
      for (int i=0; i<index_map_j.size(); i++) {
          reverse_map_j[index_map_j[i]] = i;
        }
      for (int i0=0; i0 < index0.nblock(); i0++) {
          int bs0 = index0.block_size(i0);
          int bo0 = index0.block_offset(i0);
          for (int i1=0; i1 < index1.nblock(); i1++) {
              int bs1 = index1.block_size(i1);
              int bo1 = index1.block_offset(i1);
              double maxval = 0.0;
              for (int j0=0; j0<bs0; j0++) {
                  int i_m = reverse_map_i[bo0+j0];
                  if (i_m == -1) continue;
                  for (int j1=0; j1<bs1; j1++) {
                      int j_m = reverse_map_j[bo1+j1];
                      if (j_m == -1) continue;
                      double val = fabs(m(i_m, j_m));
                      if (val > maxval) maxval = val;
                    }
                }
              if (maxval >= bound) {
                  BlockInfo<2> bi;
                  bi.block(0) = i0;
                  bi.block(1) = i1;
                  double *data = smaa.add_unallocated_block(bi).second;
                }
            }
        }
      pack_submatrix_into_array(m,smaa,index_map_i,index_map_j);
    }

    /** Packs a vector into an empty sma2::Array.  Class T must have a
        double operator()(int) member.  The nonzero blocks of the array
        will be determined. The vector does not have the full dimensions
        represented by the ranges of the array.  The index_map specifies
        how to convert the vector's indices to the array's indices. */
    template <class T>
    void pack_subvector_into_empty_array(const T &m, Array<1> &smaa,
                                         double bound,
                                         const std::vector<int> &index_map) {
      smaa.clear();
      const Range &index0 = smaa.index(0);
      std::vector<int> reverse_map(index0.nindex());
      std::fill(reverse_map.begin(), reverse_map.end(), -1);
      for (int i=0; i<index_map.size(); i++) {
          reverse_map[index_map[i]] = i;
        }
      for (int i0=0; i0 < index0.nblock(); i0++) {
          int bs0 = index0.block_size(i0);
          int bo0 = index0.block_offset(i0);
          double maxval = 0.0;
          for (int j0=0; j0<bs0; j0++) {
              int i_m = reverse_map[bo0+j0];
              if (i_m == -1) continue;
              double val = fabs(m(i_m));
              if (val > maxval) maxval = val;
            }
          if (maxval >= bound) {
              BlockInfo<1> bi;
              bi.block(0) = i0;
              double *data = smaa.add_unallocated_block(bi).second;
            }
        }
      pack_subvector_into_array(m,smaa,index_map);
    }

    /** Packs a sma2::Array into a matrix. Class T must have a double
        operator()(int,int) member. Matrix must have been allocated
        beforehand */
    template <class T>
    void pack_array_into_matrix(const Array<2> &smaa, T &m) {
      m.assign(0.0); // Initialize matrix to zero (not all blocks will be assigned below)
      const Array<2>::blockmap_t &Abm = smaa.blockmap();
      for (Array<2>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<2> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          int off1 = smaa.index(1).block_offset(Abi.block(1));
          BlockIter<2> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              m(off0+j.index(0), off1+j.index(1)) = dat[j.offset()];
            }
        }
    }

    /** Packs a sma2::Array into a vector. Class T must have a double
        operator()(int) member. The vector must have been allocated
        beforehand */
    template <class T>
    void pack_array_into_vector(const Array<1> &smaa, T &m) {
      m.assign(0.0);
      const Array<1>::blockmap_t &Abm = smaa.blockmap();
      for (Array<1>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<1> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          BlockIter<1> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              m(off0+j.index(0)) = dat[j.offset()];
            }
        }
    }

    /** Packs a sma2::Array into a subvector. Class T must have a double
        operator()(int) member. The subvector must have been allocated
        beforehand. The index map is used to convert Array indices into
        subvector indices, that is index_map[array_index_i] =
        submatrix_index_i. */
    template <class T>
    void pack_array_into_subvector(const Array<1> &smaa, T &m,
                                   const std::vector<int> &index_map) {
      m.assign(0.0);
      const Array<1>::blockmap_t &Abm = smaa.blockmap();
      for (Array<1>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<1> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          BlockIter<1> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              m(index_map[off0+j.index(0)]) = dat[j.offset()];
            }
        }
    }

    /** Packs a sma2::Array into a submatrix. Class T must have a double
        operator()(int,int) member. The submatrix must have been allocated
        beforehand.  The index maps are used to convert Array indices into
        submatrix indices, that is index_map_i[array_index_i] =
        submatrix_index_i. */
    template <class T>
    void pack_array_into_submatrix(const Array<2> &smaa, T &m,
                                   const std::vector<int> &index_map_i,
                                   const std::vector<int> &index_map_j) {
      m.assign(0.0); // Initialize matrix to zero (not all blocks will be assigned below)
      const Array<2>::blockmap_t &Abm = smaa.blockmap();
      for (Array<2>::blockmap_t::const_iterator i = Abm.begin();
           i != Abm.end();
           i++) {
          const BlockInfo<2> &Abi = i->first;
          int off0 = smaa.index(0).block_offset(Abi.block(0));
          int off1 = smaa.index(1).block_offset(Abi.block(1));
          BlockIter<2> j(smaa.indices(), Abi);
          double *dat = i->second;
          for (j.start(); j.ready(); j++) {
              m(index_map_i[off0+j.index(0)], index_map_j[off1+j.index(1)]) = dat[j.offset()];
            }
        }
    }

    /** Packs two electron integrals into array.  The array must
        already be initialized to have nonzero blocks specified.  */
    void pack_2e_integrals_into_array(Array<4> & array,
                                      const sc::Ref<sc::TwoBodyInt> &tbint);

    /** Packs two electron integrals into array.  The array must
        already be initialized to have desired blocks specified.
        Each index must be blocked by shells.  */
    void pack_2e_integrals_into_shell_blocked_array(Array<4> & array,
                                      const sc::Ref<sc::TwoBodyInt> &tbint);

}
}

#include "arraydef.h"
#include "blockinfodef.h"
#include "contract.h"
#include "contractpartdef.h"

#endif

