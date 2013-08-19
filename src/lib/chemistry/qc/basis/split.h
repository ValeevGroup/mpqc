//
// split.h
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

#ifndef _chemistry_qc_basis_split_h
#define _chemistry_qc_basis_split_h

#include <chemistry/qc/basis/gaussbas.h>

namespace sc {

/** The SplitBasisSet class is used to split a basis set's
    contractions into multiple shells.
*/
class SplitBasisSet: public GaussianBasisSet {

  protected:

    void split(const Ref<GaussianBasisSet>&,
               std::string name = std::string());

  public:

    /** The KeyVal constructor.  It does not call the basis class
        KeyVal constructor.

        <dl>

        <dt><tt>basis</tt><dd> The gives the GaussianBasisSet object
        that will be split to form this basis set. If <tt>basis</tt> is not given,
        the constructor will try to construct
        a GaussianBasisSet object from the current input and split it.


        </dl>

        */
    SplitBasisSet(const Ref<KeyVal>&);

    /** Split the given GaussianBasisSet object.
     *  @param name optional name
     *  */
    SplitBasisSet(const Ref<GaussianBasisSet>&,
                  std::string name = std::string());

    SplitBasisSet(StateIn&);

    void save_data_state(StateOut&);

    // set to 1 to debug
    static int debug() { return 0; }

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
