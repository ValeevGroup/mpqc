//
// symmint.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _chemistry_qc_integral_symmint_h
#define _chemistry_qc_integral_symmint_h

#include <util/state/state.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/petite.h>

namespace sc {

// //////////////////////////////////////////////////////////////////////////

/// Iterator over symmetry unique shell pairs
class SymmOneBodyIntIter : public OneBodyIntIter {
  protected:
    Ref<PetiteList> pl;

  public:
    SymmOneBodyIntIter(const Ref<OneBodyInt>&, const Ref<PetiteList>&);
    ~SymmOneBodyIntIter();

    void start(int ist=0, int jst=0, int ien=0, int jen=0);
    void next();

    double scale() const;

    bool cloneable() const;
    Ref<OneBodyIntIter> clone();
};

/// Iterator over symmetry unique shell quartets
class SymmTwoBodyIntIter : public TwoBodyIntIter {
  protected:
    Ref<PetiteList> pl;

  public:
    SymmTwoBodyIntIter(const Ref<TwoBodyInt>&, const Ref<PetiteList>&);
    ~SymmTwoBodyIntIter();

    void start();
    void next();

    double scale() const;
};

/// Iterator over symmetry unique shell pairs
class SymmTwoBodyTwoCenterIntIter : public TwoBodyTwoCenterIntIter {
  protected:
    Ref<PetiteList> pl;

  public:
    SymmTwoBodyTwoCenterIntIter(const Ref<TwoBodyTwoCenterInt>&, const Ref<PetiteList>&);
    ~SymmTwoBodyTwoCenterIntIter();

    void start(int ist=0, int jst=0, int ien=0, int jen=0);
    void next();

    double scale() const;

    bool cloneable() const;
    Ref<TwoBodyTwoCenterIntIter> clone();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
