//
// scfden.h
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

#ifndef _chemistry_qc_scf_scfops_h
#define _chemistry_qc_scf_scfops_h

#include <math/scmat/elemop.h>
#include <math/scmat/blocked.h>

#include <chemistry/qc/scf/scf.h>

namespace sc {

class SCFEnergy : public SCElementOp2 {
  private:
    double eelec;
    int deferred_;
    
  public:
    SCFEnergy();
    ~SCFEnergy();

    int has_collect();
    void defer_collect(int h);
    void collect(const Ref<MessageGrp>&grp);
    double result();
    void reset();

    void process(SCMatrixBlockIter&i, SCMatrixBlockIter&j);
};

class LevelShift : public BlockedSCElementOp {
  protected:
    SCF *scf_;
    double shift;

  public:
    LevelShift(SCF*);
    ~LevelShift();

    int has_side_effects();
    void set_shift(double);
    
    void process(SCMatrixBlockIter&);
};

class ALevelShift : public LevelShift {
  public:
    ALevelShift(SCF*);
    ~ALevelShift();
    void process(SCMatrixBlockIter&);
};

class BLevelShift : public LevelShift {
  public:
    BLevelShift(SCF*);
    ~BLevelShift();
    void process(SCMatrixBlockIter&);
};

// MO lagrangian
//       c  o  v
//  c  |FC|FC| 0|
//     ----------
//  o  |FC|FO| 0|
//     ----------
//  v  | 0| 0| 0|
//
class MOLagrangian : public BlockedSCElementOp2 {
  private:
    SCF *scf_;

  public:
    MOLagrangian(SCF* s);
    ~MOLagrangian();

    int has_side_effects();

    void process(SCMatrixBlockIter& bi1, SCMatrixBlockIter& bi2);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
