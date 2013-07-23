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

#ifndef _chemistry_qc_basis_dercent_h
#define _chemistry_qc_basis_dercent_h

#include <chemistry/qc/basis/basis.h>
#include <vector>

namespace sc {

/** DerivCenters keeps track the centers that
    derivatives are taken with respect to.
    DerivCenters will be used to describe a set of geometric
    derivative integral shell-sets to specify which derivatives
    are present in the set.
  */
class DerivCenters {
  private:
    std::vector< std::pair<int, int> > centers_;   // sequence of {center,atom index} pairs (center = 0 .. 3)
    std::pair<int, int> omitted_center_;           // points out which center is skipped
  public:
    /// Construct a new, empty DerivCenters object.
    DerivCenters();

    /// @name Routines to modify DerivCenters
    //@{
    /// Clear the list of centers.
    void clear();
    /** Add a center using the atom number.
     * @param center The center number (between 0 and 3 inclusive).
     * @param atom   The center within a GaussianBasisSet.
     */
    void add_center(int center, int atom);
    /** Add the omitted center using the atom number.
     * @param center The center number (between 0 and 3 inclusive).
     * @param atom   The center within a GaussianBasisSet.
     */
    void add_omitted(int center, int atom);
    //@}

    /// @name Members to query DerivCenters
    //@{
    /// The number of centers for which derivatives have been computed.
    int n() const { return centers_.size(); }
    /** @param i The computed center index (between 0 and n() - 1, inclusive).
     *  @return The center number (between 0 and 3, inclusive).
     */
    int center(int i) const { return centers_.at(i).first; }
    /** @param i The computed center index (between 0 and n() - 1, inclusive).
     *  @return The atom number.
     */
    int atom(int i) const { return centers_.at(i).second; }
    /** @return 1 if there is an omitted center, otherwise 0.
     */
    bool has_omitted_center() const { return omitted_center_.first >= 0; }
    /** @return The center for which integrals where not computed.
     */
    int omitted_center() const { return omitted_center_.first; }
    /** @return The atom that is omitted from the integral buffer.
     */
    int omitted_atom() const { return omitted_center_.second; }
    //}@
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
