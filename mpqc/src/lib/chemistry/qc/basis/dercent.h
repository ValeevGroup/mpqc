//
// dercent.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_basis_dercent_h
#define _chemistry_qc_basis_dercent_h

#include <chemistry/qc/basis/basis.h>

namespace sc {

/** DerivCenters keeps track the centers that
    derivatives are taken with respect to. */
class DerivCenters {
  private:
    int center_[4];
    int atom_[4];
    int ncenter_;
    int omitted_center_;
    int omitted_atom_;
  public:
    /** These are used by the Integral specializations to
        initializes the DerivCenters structure. */
    DerivCenters();
    /// Clear the list of centers.
    void clear();
    /// Add a center for which derivatives will be computed.
    void add_center(int center, const Ref<GaussianBasisSet> &, int shell);
    /// Add a center for which derivatives will not be computed.
    void add_omitted(int center, const Ref<GaussianBasisSet> &, int shell);
    /// Add a center for which derivatives will be computed.
    void add_center(int center, int atom);
    /// Add a center for which derivatives will not be computed.
    void add_omitted(int center, int atom);

    /// The number of unique centers minus one.
    int n() const { return ncenter_; }
    /// The center number.
    int center(int i) const { return center_[i]; }
    /// The atom number.
    int atom(int i) const { return atom_[i]; }
    /// The center that is omitted from the integral buffer.
    int omitted_center() const { return omitted_center_; }
    /// Returns 1 if there is an omitted center.
    int has_omitted_center() const { return omitted_center_ >= 0; }
    /// The atom that is omitted from the integral buffer.
    int omitted_atom() const { return omitted_atom_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
