//
// accum.h
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

#ifndef _chemistry_qc_wfn_accum_h
#define _chemistry_qc_wfn_accum_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/scmat/elemop.h>
#include <math/scmat/block.h>
#include <math/scmat/blkiter.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/integral.h>

////////////////////////////////////////////////////////////////////////////

// computes the density independent part of H
class AccumDIH: public SavableState {
#   define CLASSNAME AccumDIH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefGaussianBasisSet basis_set_;
    RefIntegral integral_;

  public:
    AccumDIH();
    AccumDIH(StateIn&);
    AccumDIH(const RefKeyVal&);
    virtual ~AccumDIH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefGaussianBasisSet&, const RefIntegral&);
    virtual void accum(const RefSymmSCMatrix& h) =0;
    virtual void done();
};
SavableState_REF_dec(AccumDIH);

////////////////////////////////////////////////////////////////////////////

// computes the density dependent part of H
class AccumDDH: public SavableState {
#   define CLASSNAME AccumDDH
#   include <util/state/stated.h>
#   include <util/class/classda.h>
  protected:
    RefGaussianBasisSet basis_set_;
    RefIntegral integral_;
    
  public:
    AccumDDH();
    AccumDDH(StateIn&);
    AccumDDH(const RefKeyVal&);
    virtual ~AccumDDH();

    void save_data_state(StateOut&);
    
    virtual void init(const RefGaussianBasisSet&, const RefIntegral&);
    virtual void accum(const RefSymmSCMatrix& h,
                       const RefSymmSCMatrix& h_open) = 0;
    virtual void done();
};
SavableState_REF_dec(AccumDDH);

class AccumNullDDH: public AccumDDH {
#   define CLASSNAME AccumNullDDH
#   define HAVE_CTOR
#   define HAVE_STATEIN_CTOR
#   define HAVE_KEYVAL_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  public:
    AccumNullDDH();
    AccumNullDDH(StateIn&);
    AccumNullDDH(const RefKeyVal&);
    ~AccumNullDDH();

    void save_data_state(StateOut&);
    
    void accum(const RefSymmSCMatrix& h, const RefSymmSCMatrix& h_open);
};

SavableState_REF_dec(AccumNullDDH);

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
