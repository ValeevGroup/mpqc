//
// pairiter.h
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _chemistry_qc_mbptr12_pairiter_h
#define _chemistry_qc_mbptr12_pairiter_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdexcept>
#include <chemistry/qc/mbptr12/moindexspace.h>

namespace sc {

/** MOPairIter gives the ordering of orbital pairs */
class MOPairIter : public RefCount {

  protected:
    bool i_eq_j_;
    int ni_;
    int nj_;
    int i_;
    int j_;
    int nij_;
    int ij_;

  public:
    /// Initialize an iterator for the given MO spaces.
    MOPairIter(const Ref<MOIndexSpace>& space_i, const Ref<MOIndexSpace>& space_j);
    virtual ~MOPairIter();

    /// Start the iteration.
    virtual void start(const int first_ij =0) =0;
    /// Move to the next pair.
    virtual void next() =0;
    /// Returns nonzero if the iterator currently hold valid data.
    virtual operator int() =0;

    /// Returns the number of functions in space i.
    int ni() { return ni_; }
    /// Returns the number of functions in space j.
    int nj() { return nj_; }
    /// Returns index i.
    int i() { return i_; }
    /// Returns index j.
    int j() { return j_; }
    /// Returns the number of pair combinations for this iterator
    int nij() { return nij_; }
    /// Returns the current iteration
    int ij() { return ij_; }
};

/** MOPairIter_SD gives the ordering of same-spin and different-spin orbital pairs.
    It iterates over all i >= j combinations (total of (ni_+1)*(ni_+2)/2 pairs). */
class MOPairIter_SD : public MOPairIter {

  int nij_aa_;
  int nij_ab_;
  int ij_aa_;
  int ij_ab_;
  int ji_ab_;

  void init_ij(const int ij) {

    if (ij<0)
      throw std::runtime_error("MOPairIter_SD::start() -- argument ij out of range");

    ij_ = 0;
    const int renorm_ij = ij%nij_;

    i_ = (int)floor((sqrt(1.0+8.0*renorm_ij) - 1.0)/2.0);
    const int i_off = i_*(i_+1)/2;
    j_ = renorm_ij - i_off;

    ij_ab_ = i_*nj_ + j_;
    ji_ab_ = j_*ni_ + i_;

    if (i_ != 0) {
      const int i_off = i_*(i_-1)/2;
      ij_aa_ = i_off + j_;
      if (i_ == j_)
        ij_aa_--;
    }
    else {
      ij_aa_ = -1;
    }
  };
  
  void inc_ij() {
    ij_++;
    if (ij_ab_ == nij_ab_-1) {
      i_ = 0;
      j_ = 0;
      ij_ab_ = 0;
      ji_ab_ = 0;
      ij_aa_ = -1;
    }
    else {
      if (i_ == j_) {
        i_++;
        j_ = 0;
        ji_ab_ = i_;
        ij_ab_ = i_*nj_;
        ij_aa_ += (i_ == j_) ? 0 : 1;
      }
      else {
        j_++;
        ji_ab_ += ni_;
        ij_ab_++;
        ij_aa_ += (i_ == j_) ? 0 : 1;
      }
    }
  };

public:
  /// Initialize an iterator for the given MO spaces.
  MOPairIter_SD(const Ref<MOIndexSpace>& space1);
  ~MOPairIter_SD();

  /// Initialize the iterator assuming that iteration will start with pair ij_offset
  void start(const int ij_offset=0)
  {
    ij_ = 0;
    init_ij(ij_offset);
  };
  
  /// Move to the next pair.
  void next() {
    inc_ij();
  };
  /// Returns nonzero if the iterator currently hold valid data.
  operator int() { return (nij_ > ij_);};

  /// Returns the number of functions in alpha-alpha space.
  int nij_aa() { return nij_aa_; }
  /// Returns the number of functions in alpha-beta space.
  int nij_ab() { return nij_ab_; }
  /** Returns compound index ij for alpha-alpha case. The i == j
      combination doesn't make sense, so ij_aa() will return -1 for such pairs. */
  int ij_aa() { return (i_ == j_) ? -1 : ij_aa_; }
  /// Returns compound index ij for alpha-beta case
  int ij_ab() { return ij_ab_; }
  /// Returns compound index ji for alpha-beta case
  int ji_ab() { return ji_ab_; }
};


/** MOPairIter_SD_neq gives the ordering of pairs of orbitals from different spaces.
It iterates over all ij combinations (total of ni_*nj_ pairs). */
class MOPairIter_SD_neq : public MOPairIter {

  int IJ_;

  void init_ij(const int ij) {

    if (ij<0)
      throw std::runtime_error("MOPairIter_SD_neq::start() -- argument ij out of range");

    IJ_ = 0;
    const int renorm_ij = ij%nij_;

    i_ = renorm_ij/nj_;
    j_ = renorm_ij - i_*nj_;

    IJ_ = i_*nj_ + j_;

  };

  void inc_ij() {
    ij_++;
    IJ_++;
    if (IJ_ == nij_) {
      i_ = 0;
      j_ = 0;
      IJ_ = 0;
    }
    else {
      if (j_ == nj_-1) {
        i_++;
        j_ = 0;
      }
      else {
        j_++;
      }
    }
  };

public:
  /// Initialize an iterator for the given MO spaces.
  MOPairIter_SD_neq(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space1);
  ~MOPairIter_SD_neq();

  /// Initialize the iterator assuming that iteration will start with pair ij_offset
  void start(const int ij_offset=0)
  {
    ij_ = 0;
    init_ij(ij_offset);
  };

  /// Move to the next pair.
  void next() {
    inc_ij();
  };
  /// Returns nonzero if the iterator currently hold valid data.
  operator int() { return (nij_ > ij_);};

  /// Returns the number of functions in alpha-alpha space.
  int nij_aa() { return nij_; }
  /// Returns the number of functions in alpha-beta space.
  int nij_ab() { return nij_; }
  /// Returns compound index ij for alpha-alpha case
  int ij_aa() { return IJ_; }
  /// Returns compound index ij for alpha-beta case
  int ij_ab() { return IJ_; }
  /// Returns compound index ij for beta-alpha case
  int ij_ba() { return IJ_; }
};


  /** This class produces MOPairIter objects */

class MOPairIterFactory {

public:
  MOPairIterFactory() {}
  ~MOPairIterFactory() {}

  /// Constructs an appropriate MOPairIter object
  Ref<MOPairIter> mopairiter(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2);
  /// Constructs an appropriate RefSCDimension object for same-spin pair
  RefSCDimension scdim_aa(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2);
  /// Constructs an appropriate RefSCDimension object for different-spin pair
  RefSCDimension scdim_ab(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2);
};

}
  
#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
