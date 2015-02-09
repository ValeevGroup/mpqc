//
// treemat.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: May 1, 2014
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
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _chemistry_qc_scf_cadf_treemat_h
#define _chemistry_qc_scf_cadf_treemat_h

#include <queue>
#include <type_traits>
#include <iterator>

#include <boost/shared_ptr.hpp>

#include <Eigen/Dense>

#include <util/misc/assert.h>
#include <util/misc/scexception.h>

#include "treemat_fwd.h"
#include "cadfclhf.h"

namespace sc { namespace cadf {

template<typename NormContainer, typename Index, typename NormValue>
class TreeBlock {
  public:

    typedef NormContainer norm_container_t;
    typedef Index index_t;
    typedef NormValue norm_value_t;

  private:

    norm_container_t norms_;
    index_t begin_index_;
    index_t end_index_;
    std::vector<boost::shared_ptr<TreeBlock>> children_;
    int atom_block_index_ = -1;

  public:

    TreeBlock() = default;

    template<typename Derived>
    TreeBlock(
        const Eigen::MatrixBase<Derived>& norms,
        index_t begin_idx, index_t end_idx
    ) : norms_(norms),
        begin_index_(begin_idx),
        end_index_(end_idx)
    { }

    template<typename Iterator>
    static boost::shared_ptr<TreeBlock>
    merge_blocks(Iterator begin, Iterator end, Index end_index, int atom_block_index = -1)
    {
      boost::shared_ptr<TreeBlock> rv(boost::make_shared<TreeBlock>());
      const Index norms_size = (*begin)->norms_.size();
      rv->norms_.resize(norms_size);
      rv->norms_ = NormContainer::Zero(norms_size);
      rv->begin_index_ = (*begin)->begin_index();
      rv->end_index_ = end_index;
      rv->atom_block_index_ = atom_block_index;

      // We can only reserve the size of the children vector ahead of time if we can get the
      //   difference between the iterators
      if(std::is_base_of<std::random_access_iterator_tag,
          typename std::iterator_traits<Iterator>::iterator_category>::value
      )
      {
        typename std::iterator_traits<Iterator>::difference_type dist = end - begin;
        rv->children_.reserve(dist);
      }

      for(auto it = begin; it != end; ++it) {
        rv->norms_.array() += (*it)->norms_.array().square();
        rv->children_.push_back(*it);
      }
      rv->norms_ = rv->norms_.cwiseSqrt();
      return rv;
    }

    index_t begin_index() const { return begin_index_; }

    index_t end_index() const { return end_index_; }

    bool is_atom_block() const { return atom_block_index_ != -1; }

    int atom_block_index() const { return atom_block_index_; }

    template <typename NormIndex>
    const norm_value_t
    norm(const NormIndex idx) const
    {
      return norms_[idx];
    }

    size_t n_children() const { return children_.size(); }
    bool is_leaf() const { return n_children() == 0; }

    const std::vector<boost::shared_ptr<TreeBlock>>& children() const { return children_; }

};

template<typename BlockType>
class TreeMatrix
{
  public:

    typedef BlockType block_t;
    typedef typename block_t::index_t index_t;
    typedef boost::shared_ptr<block_t> block_ptr_t;
    typedef std::deque<block_ptr_t> block_container_t;

  private:

    // Blocks ordered from root to last leaf
    block_container_t blocks_;

  public:

    const block_container_t& blocks() const { return blocks_; }

    typename block_container_t::const_iterator
    breadth_first_begin() const {
      return blocks_.cbegin();
    }

    typename block_container_t::const_iterator
    breadth_first_end() const {
      return blocks_.cend();
    }

    const block_t root() const { return blocks_.front(); }

    template<typename Derived>
    TreeMatrix(
        const Eigen::MatrixBase<Derived>& m,
        GaussianBasisSet* basis,
        std::vector<int> block_requirements = { SameAngularMomentum|SameCenter, SameCenter },
        int max_children = 3
    )
    {
      // These use std::deque objects since we need the iterators, but they act like stacks
      std::deque<block_ptr_t> curr_blocks;
      std::deque<block_ptr_t> new_blocks;
      for(const auto&& ish : shell_range(basis))
      {
        auto blk = boost::make_shared<block_t>(m.row(ish), ish.bfoff, ish.bfoff+ish.nbf);
        new_blocks.push_front(blk);
      }

      for(auto req : block_requirements) {
        curr_blocks.clear();
        while(!new_blocks.empty()) {
          curr_blocks.push_front(new_blocks.front());
          blocks_.push_front(new_blocks.front());
          new_blocks.pop_front();
        }

        auto blk_iter = curr_blocks.cbegin();
        for(const auto&& iblk : shell_block_range(basis, 0, 0, NoLastIndex, req, NoMaximumBlockSize)) {
          auto blk_start = blk_iter;
          index_t end_index = (*blk_start)->end_index();
          while(blk_iter != curr_blocks.cend()
              and (*blk_iter)->begin_index() < iblk.bfoff + iblk.nbf
          ) {
            end_index = (*blk_iter)->end_index();
            ++blk_iter;
          }
          new_blocks.push_front(block_t::merge_blocks(blk_start, blk_iter, end_index,
              req==SameCenter ? iblk.center : -1
          ));
        }
      }

      while(new_blocks.size() > 1) {

        curr_blocks.clear();
        while(!new_blocks.empty()) {
          curr_blocks.push_front(new_blocks.front());
          blocks_.push_front(new_blocks.front());
          new_blocks.pop_front();
        }

        auto blk_iter = curr_blocks.cbegin();
        while(blk_iter != curr_blocks.cend()) {
          auto blk_start = blk_iter;
          index_t end_index = (*blk_start)->end_index();
          for(int ichild = 0; ichild < max_children and blk_iter != curr_blocks.end(); ++ichild) {
            end_index = (*blk_iter)->end_index();
            ++blk_iter;
          }
          new_blocks.push_front(block_t::merge_blocks(blk_start, blk_iter, end_index));
        }

      }

      blocks_.push_front(new_blocks.front());

    }


};

template<
  typename Index=uli,
  typename LeftBlockPtr=typename TreeMatrix<>::block_ptr_t,
  typename RightBlockPtr=typename TreeMatrix<>::block_ptr_t
>
class ProductBlock {

  private:

    double norm_;
    Index begin_index_;
    Index end_index_;
    LeftBlockPtr left_block_;
    RightBlockPtr right_block_;

  public:

    template<typename LeftIndex, typename RightIndex>
    ProductBlock(
        const LeftBlockPtr& left, LeftIndex left_index,
        const RightBlockPtr& right, RightIndex right_index
    ) : norm_(left->norm(left_index) * right->norm(right_index)),
        begin_index_(left->begin_index()),
        end_index_(left->end_index()),
        left_block_(left), right_block_(right)
    {
      MPQC_ASSERT(left->begin_index() == right->begin_index());
      MPQC_ASSERT(left->end_index() == right->end_index());
      MPQC_ASSERT(left->n_children() == right->n_children());
    }

    const Index size() const { return end_index_ - begin_index_; }
    const Index begin_index() const { return begin_index_; }
    const Index end_index() const { return end_index_; }
    const double norm() const { return norm_; }

    bool operator<(const ProductBlock& other) const {
      return norm_ / size() < other.norm_ / other.size();
    }

    template<typename LeftIndex, typename RightIndex>
    void
    enqueue_children(
        std::queue<ProductBlock>& queue,
        const LeftIndex left_index,
        const RightIndex right_index
    ) const
    {
      auto left_iter = left_block_->children().cbegin();
      auto right_iter = right_block_->children().cbegin();
      for(; left_iter != left_block_->children().end(); ++left_iter, ++right_iter) {
        queue.emplace(*left_iter, left_index, *right_iter, right_index);
      }
    }

    bool is_leaf() const {
      return left_block_->is_leaf();
    }

    bool is_atom_block() const {
      return left_block_->is_atom_block();
    }

    int atom_block_index() const {
      return left_block_->atom_block_index();
    }


};

namespace detail {

  template<typename T>
  struct begin_index_less
  {
      bool operator()(const T& a, const T& b) const {
        return a.begin_index() < b.begin_index();
      }
  };

}

template<typename LeftTreeType, typename RightTreeType, typename LeftIndex, typename RightIndex>
inline std::vector<std::pair<
    typename LeftTreeType::index_t,
    typename RightTreeType::index_t
>>
relevant_product_ranges(
    const LeftTreeType& left,  LeftIndex left_index,
    const RightTreeType& right, RightIndex right_index,
    double thresh, int exclude_center = -1,
    // Note: Using this option is possibly more stable
    //   but could result in a *slightly* stronger dependence on
    //   ordering of atoms.
    bool guarantee_min_error=false
)
{
  typedef ProductBlock<
      typename LeftTreeType::index_t,
      typename LeftTreeType::block_ptr_t,
      typename RightTreeType::block_ptr_t
    > product_block_t;
  typedef std::vector<std::pair<
      typename LeftTreeType::index_t,
      typename RightTreeType::index_t
    >> return_t;

  std::priority_queue<product_block_t> discarded_blocks;
  std::set<product_block_t, detail::begin_index_less<product_block_t>> sig_blocks;
  std::queue<product_block_t> working_blocks;
  double curr_error_squared = 0.0;

  auto left_iter = left.breadth_first_begin();
  auto right_iter = right.breadth_first_begin();

  working_blocks.emplace(*left_iter, left_index, *right_iter, right_index);

  while(!working_blocks.empty()) {
    // TODO think about how to avoid swapping left and right block lists in and out of cache (perhaps by prefetching?)
    const auto& curr_block = working_blocks.front();

    if(exclude_center != -1 and curr_block.is_atom_block() and curr_block.atom_block_index() == exclude_center) {
      // Discard the block; don't even put it in the discarded blocks
    }
    else if(curr_block.norm() < thresh) {
      if(guarantee_min_error) {
        discarded_blocks.push(curr_block);
      }
      curr_error_squared += curr_block.norm() * curr_block.norm();
    }
    else {
      if(curr_block.is_leaf()) {
        sig_blocks.insert(curr_block);
      }
      else {
        curr_block.enqueue_children(working_blocks, left_index, right_index);
      }
    }
    working_blocks.pop();
  }

  if(guarantee_min_error) {
    throw FeatureNotImplemented("guarantee_min_error", __FILE__, __LINE__);
  }

  // Form the largest contiguous blocks that we can
  return_t rv;
  for(const auto& sig_block : sig_blocks) {
    if(rv.empty() or rv.back().second != sig_block.begin_index()) {
      rv.emplace_back(sig_block.begin_index(), sig_block.end_index());
    }
    else {
      rv.back().second = sig_block.end_index();
    }
  }

  return rv;
}


}} // end namespaces

#endif /* _chemistry_qc_scf_cadf_treemat_h */
