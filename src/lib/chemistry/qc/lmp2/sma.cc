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

#include <algorithm>

#include <chemistry/qc/lmp2/sma.h>

namespace sc {

namespace sma2 {

    Data::Data():
      ndata_(0)
    {
      default_chunksize_ = 4096;
    }

    Data::~Data()
    {
      for (memmap_t::iterator iter = memmap_.begin();
           iter != memmap_.end(); iter++) {
          delete[] iter->second.first;
        }
    }

    double *
    Data::allocate(long size)
    {
      double *data = new double[size];
      memitermap_[data]
          = memmap_.insert(std::make_pair(0,std::make_pair(data,size)));
      ndata_ += size;
      return data;
    }

    void
    Data::deallocate(double *dat)
    {
      memitermap_t::iterator memiter = memitermap_.find(dat);
      if (memiter == memitermap_.end()) {
          throw std::runtime_error("Data::deallocate: data not found");
        }
      ndata_ -= memiter->second->second.second;
      memmap_.erase(memiter->second);
      memitermap_.erase(memiter);
      delete[] dat;
    }

    double *
    Data::data() const
    {
      if (memmap_.size() == 1) {
          return memmap_.begin()->second.first;
        }
      throw std::runtime_error("Data::data not implemented");
    }

    //////////////////////////////////////////////////////////

    IndexList::IndexList(): indices_(0)
    {
    }

    IndexList::IndexList(int a): indices_(1)
    {
      indices_[0] = a;
    }

    IndexList::IndexList(int a,int b): indices_(2)
    {
      indices_[0] = a;
      indices_[1] = b;
    }

    IndexList::IndexList(int a,int b,int c): indices_(3)
    {
      indices_[0] = a;
      indices_[1] = b;
      indices_[2] = c;
    }

    IndexList::IndexList(int a,int b,int c,int d): indices_(4)
    {
      indices_[0] = a;
      indices_[1] = b;
      indices_[2] = c;
      indices_[3] = d;
    }

    IndexList::IndexList(const std::vector<int> &v): indices_(v)
    {
    }

    IndexList::IndexList(const IndexList &i):
      indices_(i.indices_)
    {
    }

    IndexList::IndexList(const IndexList &i1, const IndexList &i2):
      indices_(i1.indices_)
    {
      append_additional_indices(i2);
    }

    void
    IndexList::append_additional_indices(const IndexList &il)
    {
      for (int i=0; i<il.n(); i++) indices_.push_back(il.i(i));
    }

    IndexList
    IndexList::reverse_mapping() const
    {
      IndexList r;
      r.indices_.resize(n());
      for (int i=0; i<n(); i++) r.indices_[indices_[i]] = i;
      return r;
    }

    IndexList
    IndexList::identity(int n)
    {
      IndexList r;
      r.indices_.resize(n);
      for (int i=0; i<n; i++) r.indices_[i] = i;
      return r;
    }

    bool
    IndexList::is_identity() const
    {
      for (int i=0; i<n(); i++) {
          if (indices_[i] != i) return false;
        }
      return true;
    }

    bool
    IndexList::is_identity_permutation() const
    {
      std::vector<int> indices_copy = indices_;
      std::sort(indices_copy.begin(), indices_copy.end());
      for (int i=0; i<indices_copy.size(); i++)
          if (indices_copy[i] != i) return false;
      return true;
    }

    void
    IndexList::print(std::ostream &o) const {
      o << "{";
      for (int i=0; i<indices_.size(); i++) {
          if (i>0) o << " ";
          o << indices_[i];
        }
      o << "}";
    }

    std::ostream &
    operator << (std::ostream &o, const IndexList &l)
    {
      l.print(o);
      return o;
    }

    //////////////////////////////////////////////////////////

    void
    Range::init_offsets()
    {
      block_offset_.clear();
      int offset = 0;
      for (int i=0; i<block_size_.size(); i++) {
          block_offset_.push_back(offset);
          offset += block_size_[i];
        }
      nindex_ = offset;

      index_to_block_.resize(nindex_);
      int index = 0;
      for (int i=0; i<block_size_.size(); i++) {
          for (int j=0; j<block_size_[i]; j++, index++) {
              index_to_block_[index] = i;
            }
        }
      range_order_.resize(function_order_.size());
      for (int i=0; i<function_order_.size(); i++) {
          range_order_[function_order_[i]] = i;
        }

      max_block_size_ = 0;
      for (int i=0; i<block_size_.size(); i++) {
          if (block_size_[i] > max_block_size_)
              max_block_size_ = block_size_[i];
        }
    }

    void
    Range::init_extent_blocking(const sc::Ref<sc::GaussianBasisSet> &bs)
    {
      function_order_.resize(bs->nbasis());
      int ifunc = 0;
      for (int i=0; i<bs->ncenter(); i++) {
          // sort the shells on this center from smallest to largest
          // spatial extent by placing extent, shell number pairs into a
          // multimap
          std::multimap<double,int> shellmap;
          for (int j=0; j<bs->nshell_on_center(i); j++) {
              shellmap.insert(
                  std::pair<double,int>(
                      bs->shell(bs->shell_on_center(i,j)).extent(1.0e-4),
                      j
                      )
                  );
            }
          // the maximum extent for the current block
          // this value is chosen to put hydrogen AO's
          // in the same block for 3-21G
          double maxextent = 6.5;
          double delta_maxextent = 4.0;
          int currentblocksize = 0;
          int width = 0;
          for (std::multimap<double,int>::iterator j = shellmap.begin();
               j != shellmap.end();
               j++) {
              int jsh = j->second;
              int ijsh = bs->shell_on_center(i,jsh);
              if (j->first > maxextent) {
                  if (currentblocksize) {
                      block_size_.push_back(currentblocksize);
                      currentblocksize = 0;
                    }
                  do {
                      maxextent += delta_maxextent;
                    } while (j->first > maxextent);
                }
              // add shell to current block
              int nfunction = bs->shell(ijsh).nfunction();
              currentblocksize += nfunction;
              int function_offset = bs->shell_to_function(ijsh);
              for (int k=0; k<nfunction; k++,ifunc++) {
                  function_order_[ifunc] = function_offset + k;
                }
            }
          if (currentblocksize) {
              block_size_.push_back(currentblocksize);
            }
        }
    }

    void
    Range::init_order(int nbasis)
    {
      function_order_.resize(nbasis);
      for (int i=0; i<nbasis; i++) function_order_[i] = i;
    }

    Range::Range(const sc::Ref<sc::GaussianBasisSet> &bs,
                 BlockingMethod bm, int blocksize)
    {
      init(bs,bm,blocksize);
    }

    Range::Range(const std::vector<int> &block_size)
    {
      nindex_ = 0;
      for (int i=0; i<block_size.size(); i++) nindex_ += block_size[i];
      block_size_ = block_size;
      init_order(nindex_);
      init_offsets();
    }

    void
    Range::init(const sc::Ref<sc::GaussianBasisSet> &bs,
                 BlockingMethod bm, int blocksize)
    {
      clear();
      if (bm == FunctionBlocking) {
          int nbasis_remaining = bs->nbasis();
          while (nbasis_remaining > 0) {
              int current_blocksize;
              if (nbasis_remaining > blocksize) current_blocksize = blocksize;
              else current_blocksize = nbasis_remaining;
              block_size_.push_back(current_blocksize);
              nbasis_remaining -= current_blocksize;
            }
          init_order(bs->nbasis());
        }
      else if (bm == ShellBlocking) {
          for (int i=0; i<bs->nshell(); i++) {
              block_size_.push_back(bs->shell(i).nfunction());
            }
          init_order(bs->nbasis());
        }
      else if (bm == AtomBlocking) {
          for (int i=0; i<bs->ncenter(); i++) {
              block_size_.push_back(bs->nbasis_on_center(i));
            }
          init_order(bs->nbasis());
        }
      else if (bm == ExtentBlocking) {
          init_extent_blocking(bs);
        }
      else {
          throw std::invalid_argument("unknown blocking scheme");
        }
      init_offsets();
    }

    Range::Range(int nindex, int block_size)
    {
      init(nindex, block_size);
    }

    void
    Range::init(int nindex, int block_size)
    {
      clear();
      int nremain = nindex;
      while (nremain > 0) {
          if (nremain > block_size) {
              block_size_.push_back(block_size);
              nremain -= block_size;
            }
          else {
              block_size_.push_back(nremain);
              nremain = 0;
            }
        }
      init_order(nindex);
      init_offsets();
    }

    Range::Range()
    {
      clear();
    }

    void
    Range::clear()
    {
      nindex_ = 0;
      block_size_.clear();
      block_offset_.clear();
      index_to_block_.clear();
      function_order_.clear();
      range_order_.clear();
    }

    bool
    Range::operator == (const Range &r) const
    {
      if (nindex_ != r.nindex_) return false;
      int nblock = this->nblock();
      if (nblock != r.nblock()) return false;
      for (int i=0; i<nblock; i++) {
          if (block_size_[i] != r.block_size_[i]) return false;
        }
      return true;
    }

    bool
    Range::all_size_one() const
    {
      for (int i=0; i<nblock(); i++) {
          if (block_size_[i] != 1) return false;
        }
      return true;
    }

    int
    Range::basis_index_to_range_index(int i) const
    {
      return function_order_[i];
    }

    int
    Range::range_index_to_basis_index(int i) const
    {
      return range_order_[i];
    }

    void
    Range::print(std::ostream&o) const
    {
      o << sc::indent
        << "sma::Range: nblock = " << nblock()
        << " nindex = " << nindex()
        << std::endl;
      o << sc::incindent;
      for (int ij=0, i=0; i<nblock(); i++) {
          o << sc::indent << sc::scprintf("%3d:",i);
          int linelength = 4;
          for (int j=0; j<block_size(i); j++,ij++) {
              if (linelength > 70) {
                  linelength = 0;
                  o << std::endl << sc::indent << "    ";
                  linelength = 4;
                }
              else linelength += 4;
              o << sc::scprintf(" %3d", range_index_to_basis_index(ij));
            }
          o << std::endl;
        }
      o << sc::decindent;
    }

    void
    Range::write(sc::StateOut&so) const {
      so.put(int(block_size_.size()));
      so.put(int(function_order_.size()));
      for (int i=0; i<block_size_.size(); i++) {
          so.put(block_size_[i]);
        }
      for (int i=0; i<function_order_.size(); i++) {
          so.put(function_order_[i]);
        }
    }

    void
    Range::read(sc::StateIn&si) {
      int nblock;
      si.get(nblock);
      block_size_.resize(nblock);
      int nbasis;
      si.get(nbasis);
      function_order_.resize(nbasis);
      for (int i=0; i<nblock; i++) {
          si.get(block_size_[i]);
        }
      for (int i=0; i<nbasis; i++) {
          si.get(function_order_[i]);
        }
      init_offsets();
    }

    //////////////////////////////////////////////////////////

}

}

