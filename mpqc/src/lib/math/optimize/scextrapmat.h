//
// scextrapmat.h
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

#ifndef _math_optimize_scextrapmat_h
#define _math_optimize_scextrapmat_h

#include <math/optimize/scextrap.h>
#include <math/scmat/matrix.h>

class SymmSCMatrixSCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrixSCExtrapData
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m;
  public:
    SymmSCMatrixSCExtrapData(StateIn&);
    SymmSCMatrixSCExtrapData(const RefSymmSCMatrix&);

    void save_data_state(StateOut&);
    
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrix2SCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrix2SCExtrapData
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m1;
    RefSymmSCMatrix m2;
  public:
    SymmSCMatrix2SCExtrapData(StateIn&);
    SymmSCMatrix2SCExtrapData(const RefSymmSCMatrix&, const RefSymmSCMatrix&);

    void save_data_state(StateOut&);
    
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrixNSCExtrapData: public SCExtrapData {
#   define CLASSNAME SymmSCMatrixNSCExtrapData
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int n_;
    RefSymmSCMatrix *m;
  public:
    SymmSCMatrixNSCExtrapData(StateIn&);
    SymmSCMatrixNSCExtrapData(int n, RefSymmSCMatrix*);

    void save_data_state(StateOut&);
    
    SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const RefSCExtrapData&);
};

class SymmSCMatrixSCExtrapError: public SCExtrapError {
#   define CLASSNAME SymmSCMatrixSCExtrapError
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    RefSymmSCMatrix m;
  public:
    SymmSCMatrixSCExtrapError(StateIn&);
    SymmSCMatrixSCExtrapError(const RefSymmSCMatrix&);

    void save_data_state(StateOut&);
    
    double error();
    double scalar_product(const RefSCExtrapError&);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
