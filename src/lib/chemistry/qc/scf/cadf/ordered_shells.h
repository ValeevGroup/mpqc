//
// ordered_shells.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: May 5, 2014
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

#ifndef _chemistry_qc_scf_cadf_ordered_shells_h
#define _chemistry_qc_scf_cadf_ordered_shells_h

#include <array>
#include <mutex>

#include <Eigen/Dense>

#include "iters.h"

namespace sc {

class OrderedShellList {

  public:

    typedef std::vector<ShellIndexWithValue> index_list;
    typedef std::unordered_set<
        ShellIndexWithValue,
        detail::hash_<ShellIndexWithValue>,
        detail::index_equal_
    > index_set;

    typedef index_list::const_iterator index_iterator;
    typedef basis_element_with_value_iterator<ShellDataWithValue, index_iterator> iterator;

  private:

    // TODO Make this a shared mutex
    mutable std::mutex insert_mtx_;
    mutable std::mutex aux_vector_mtx_;

    index_list indices_;
    index_set idx_set_;
    bool sorted_ = false;
    bool sort_by_value_ = true;

    Eigen::VectorXd aux_vector_;
    bool aux_vector_initialized_ = false;
    bool aux_set_sorted_ = false;
    double aux_value_ = 0.0;

    std::mutex aux_set_mtx_;
    std::unordered_set<int> aux_set_;
    std::vector<int> aux_indices_;

  public:

    GaussianBasisSet* basis_ = 0;
    GaussianBasisSet* dfbasis_ = 0;

    int nbf = 0;

    // An auxiliary value to tag along with the list

    explicit OrderedShellList(bool sort_by_value = true)
      : sort_by_value_(sort_by_value)
    {
      aux_value_ = 0.0;
      aux_vector_initialized_ = false;
    }

    template<typename Iterable>
    OrderedShellList(
        const Iterable& indices_in,
        GaussianBasisSet* basis,
        GaussianBasisSet* dfbasis = 0
    )
      : sort_by_value_(false),
        basis_(basis), dfbasis_(dfbasis)
    {
      std::copy(indices_in.begin(), indices_in.end(), std::back_inserter(indices_));
      sort(false);
      aux_value_ = 0.0;
    }

    OrderedShellList(const OrderedShellList& other)
      : indices_(other.indices_),
        idx_set_(other.idx_set_),
        sorted_(other.sorted_),
        sort_by_value_(other.sort_by_value_),
        nbf(other.nbf)
    {
      aux_value_ = other.aux_value_;
      aux_vector_initialized_ = other.aux_vector_initialized_;
    }

    void add_to_aux_value(const double add_val) {
      std::lock_guard<std::mutex> lg(aux_vector_mtx_);
      aux_value_ += add_val;
    }

    template <typename Derived>
    void add_to_aux_value_vector(const double add_val, const Eigen::MatrixBase<Derived>& to_add) {
      std::lock_guard<std::mutex> lg(aux_vector_mtx_);
      if(not aux_vector_initialized_) {
        aux_vector_.resize(to_add.rows());
        aux_vector_ = decltype(aux_vector_)::Zero(to_add.rows());
        aux_vector_initialized_ = true;
      }
      aux_value_ += add_val;
      aux_vector_.noalias() += to_add;
    }

    void set_aux_value(const double add_val) {
      aux_value_ = add_val;
    }

    void insert_in_aux_set(int item) {
      std::lock_guard<std::mutex> lg(aux_set_mtx_);
      aux_set_.insert(item);
      aux_set_sorted_ = false;
    }

    //// Note: not thread safe, do not call if insertion may happen simultaneously
    const std::vector<int>& get_aux_indices() const {
      assert(aux_set_sorted_);
      return aux_indices_;
    }

    //// Note: not thread safe, do not call if insertion may happen simultaneously
    const std::unordered_set<int>& get_aux_set() const {
      return aux_set_;
    }

    void sort_aux_set() {
      std::lock_guard<std::mutex> lg(aux_set_mtx_);
      assert(aux_indices_.empty());
      std::copy(aux_set_.begin(), aux_set_.end(), std::back_inserter(aux_indices_));
      std::sort(aux_indices_.begin(), aux_indices_.end());
      aux_set_sorted_ = true;
    }

    //// Note: not thread safe, do not call if insertion may happen simultaneously
    bool aux_set_contains(int item) const {
      return aux_set_.find(item) != aux_set_.end();
    }

    double get_aux_value() const {
      return aux_value_;
    }

    const decltype(aux_vector_)& get_aux_vector() const {
      assert(aux_vector_initialized_);
      return aux_vector_;
    }

    void set_basis(GaussianBasisSet* basis, GaussianBasisSet* dfbasis = 0)
    {
      basis_ = basis;
      if(dfbasis) dfbasis_ = dfbasis;
    }

    void insert(const ShellData& ish, double value = 0, int nbf_in = 0) {
      //----------------------------------------//
      assert(basis_ == 0 || ish.basis == basis_);
      if(basis_ == 0) basis_ = ish.basis;
      assert(dfbasis_ == 0 || ish.dfbasis == dfbasis_);
      if(dfbasis_ == 0) dfbasis_ = ish.dfbasis;
      //----------------------------------------//
      std::lock_guard<std::mutex> lg(insert_mtx_);
      sorted_ = false;
      ShellIndexWithValue insert_val(ish, value);
      const auto& found = idx_set_.find(insert_val);
      if(found != idx_set_.end()) {
        const double old_val = found->value;
        if(value > old_val) found->value = value;
      }
      else {
        idx_set_.insert(insert_val);
        nbf += nbf_in;
      }
    }

    size_t size() const {
      if(sorted_) return indices_.size();
      else return idx_set_.size();
    }

    void set_sort_by_value(bool new_srt=true) {
      if(new_srt != sort_by_value_) {
        sort_by_value_ = new_srt;
        sorted_ = false;
      }
    }

    double value_for_index(int index) {
      if(idx_set_.size() < indices_.size()) {
        idx_set_.clear();
        for(auto&& idx : indices_) {
          idx_set_.insert(idx);
        }
      }
      std::lock_guard<std::mutex> lg(insert_mtx_);
      ShellIndexWithValue find_val(index);
      const auto& found = idx_set_.find(find_val);
      if(found != idx_set_.end()) {
        return found->value;
      }
      else {
        throw AlgorithmException("Index not found in list", __FILE__, __LINE__);
      }
    }

    void sort(bool transfer_idx_set = true) {
      std::lock_guard<std::mutex> lg(insert_mtx_);
      if(transfer_idx_set) {
        assert(idx_set_.size() >= indices_.size());
        indices_.clear();
        std::copy(idx_set_.begin(), idx_set_.end(), std::back_inserter(indices_));
      }
      if(sort_by_value_) {
        std::sort(indices_.begin(), indices_.end(),
            [](const ShellIndexWithValue& a, const ShellIndexWithValue& b){
              return a.value > b.value or (a.value == b.value and a.index < b.index);
            }
        );
      }
      else {
        std::sort(indices_.begin(), indices_.end(),
            [](const ShellIndexWithValue& a, const ShellIndexWithValue& b){
              return a.index < b.index;
            }
        );
      }
      sorted_ = true;
    }

    void acquire_and_sort(
        ShellIndexWithValue* new_data,
        size_t n_new_data
    )
    {
      assert(size() == 0);
      std::copy(new_data, new_data + n_new_data, std::back_inserter(indices_));
      sort(false);
    }

    template <typename ValueIterable>
    void acquire_and_sort(
        const ValueIterable& val_iter,
        double cutoff = 0.0
    )
    {
      // There is probably a faster way to do this, particularly for Eigen data structures
      int index = 0;
      indices_.clear();
      idx_set_.clear();
      sorted_ = false;
      for(auto&& val : val_iter) {
        if(val > cutoff) {
          indices_.emplace_back(index++, val);
        }
      }
      sort(false);
    }

    void acquire_and_sort(
        const double* vals,
        size_t n_vals,
        double cutoff = 0.0,
        bool sort_by_val = true
    )
    {
      indices_.clear();
      idx_set_.clear();
      sorted_ = false;
      for(size_t idx = 0; idx < n_vals; ++idx) {
        if(vals[idx] > cutoff) {
          indices_.emplace_back(idx, vals[idx]);
        }
      }
      set_sort_by_value(sort_by_val);
      sort(false);
    }

    template<typename IndexIterable>
    void acquire_and_sort(
        const IndexIterable& idxs_in,
        const double* vals_in,
        double cutoff = 0.0,
        bool sort_by_val = true
    )
    {
      indices_.clear();
      idx_set_.clear();
      typename IndexIterable::iterator idxiter = idxs_in.begin();
      int val_off = 0;
      for(; idxiter != idxs_in.end(); ++idxiter, ++val_off) {
        if(vals_in[val_off] > cutoff) {
          indices_.emplace_back(*idxiter, vals_in[val_off]);
        }
      }
      set_sort_by_value(sort_by_val);
      sort(false);
    }

    template <typename IndexIterable>
    OrderedShellList
    intersection_with(const IndexIterable& other_indices) {
      if(!sorted_) {
        sort();
      }
      OrderedShellList rv(false);
      std::copy(indices_.begin(), indices_.end(), std::back_inserter(rv.indices_));
      if(sort_by_value_) {
        rv.sort(false);
      }
      IndexIterable other_copy;
      std::copy(other_indices.begin(), other_indices.end(), std::back_inserter(other_copy));
      std::sort(other_copy.begin(), other_copy.end(), std::less<int>());

      // do the intersection
      std::vector<ShellIndexWithValue> intersect;
      std::set_intersection(
          rv.indices_.begin(), rv.indices_.end(),
          other_copy.begin(), other_copy.end(),
          std::back_inserter(intersect), std::less<int>()
      );

      rv.indices_.clear();
      std::move(intersect.begin(), intersect.end(), std::back_inserter(rv.indices_));

      if(sort_by_value_) {
        rv.sort_by_value_ = true;
        rv.sorted_ = false;
        rv.sort(false);
      }

      rv.basis_ = basis_;
      rv.dfbasis_ = dfbasis_;

      return rv;

    }

    iterator begin() const
    {
      assert(sorted_);
      return iterator(
          basis_, dfbasis_,
          index_begin()
      );
    }

    index_iterator index_begin() const
    {
      assert(sorted_);
      return indices_.cbegin();
    }

    iterator end() const
    {
      assert(sorted_);
      return iterator(
          basis_, dfbasis_,
          index_end()
      );
    }

    index_iterator index_end() const
    {
      assert(sorted_);
      return indices_.cend();
    }

    const index_list& unsorted_indices() {
      std::lock_guard<std::mutex> lg(insert_mtx_);
      assert(!sorted_);
      if(indices_.size() != idx_set_.size()) {
        indices_.clear();
        std::copy(idx_set_.begin(), idx_set_.end(), std::back_inserter(indices_));
      }
      return indices_;
    }

    void clear()
    {
      indices_.clear();
      idx_set_.clear();
    }

    friend std::ostream&
    operator <<(std::ostream& out, const OrderedShellList& list);

};

////////////////////////////////////////////////////////////////////////////////

typedef enum {
  L3 = 0,
  L3_star = 1,
  LB = 2,
  Ld_over = 3,
  Ld_under = 4
} LinKListName;

class LinKListGroup
{
    std::array<OrderedShellList, 5> lists_;

  public:

    OrderedShellList& operator[](LinKListName name) { return lists_[name]; }

};

////////////////////////////////////////////////////////////////////////////////

inline range_of_shell_blocks<typename OrderedShellList::index_iterator>
shell_block_range(
    const OrderedShellList& shlist,
    int requirements=SameCenter,
    int target_size=DEFAULT_TARGET_BLOCK_SIZE
)
{
  return boost::make_iterator_range(
      shell_block_iterator<typename OrderedShellList::index_iterator>(
          shlist.index_begin(), shlist.index_end(),
          shlist.basis_, shlist.dfbasis_, requirements, target_size
      ),
      shell_block_iterator<typename OrderedShellList::index_iterator>(
          shlist.index_end(), shlist.index_end(),
          shlist.basis_, shlist.dfbasis_, requirements, target_size
      )
  );
}

////////////////////////////////////////////////////////////////////////////////


} // end namespace sc

#endif /* _chemistry_qc_scf_cadf_ordered_shells_h */
