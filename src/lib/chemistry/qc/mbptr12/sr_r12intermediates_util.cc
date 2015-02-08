//
// sr_r12intermediates_util.cc
//
// Copyright (C) 2013 Edward Valeev
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
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma implementation
#endif

#include <mpqc_config.h>

#if defined(MPQC_NEW_FEATURES)

#include <chemistry/qc/mbptr12/sr_r12intermediates.h>

namespace sc {

  template <>
  DA4_Tile<double>::operator TA::Tensor<double>() const {

    eval_type tile(owner_->trange().make_tile_range(index_));

    size_t size34 = tile.range().size()[2] * tile.range().size()[3];

    for(size_t e1=tile.range().start()[0], e12=0; e1!=tile.range().finish()[0]; ++e1) {
      for(size_t e2=tile.range().start()[1]; e2!=tile.range().finish()[1]; ++e2, ++e12) {

        darray4_->retrieve_pair_subblock(e1, e2, te_type_,
            tile.range().start()[2], tile.range().finish()[2],
            tile.range().start()[3], tile.range().finish()[3],
            tile.data() + e12*size34);

      }
    }

    return tile;
  }

  template <>
  DA4_Tile34<double>::operator TA::Tensor<TA::Tensor<double> >() const {

    // make ordinary ranges needed to make 34 tiles
    std::vector<size_t> i3i4_start(2, 0);
    std::vector<size_t> i3i4_finish(2); i3i4_finish[0] = darray4_->nx(); i3i4_finish[1] = darray4_->ny();
    TA::Range i3i4_range(i3i4_start, i3i4_finish);

    eval_type tile(owner_->trange().make_tile_range(index_), eval_type::value_type(i3i4_range));

    darray4_->retrieve_pair_block(index_[0], index_[1], te_type_,
        tile.data()->data());
    darray4_->release_pair_block(index_[0], index_[1], te_type_);

    return tile;
  }

}; // namespace sc

#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
