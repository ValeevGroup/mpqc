//
// effh.h --- definition of the effective fock builder classes
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

#ifndef _chemistry_qc_scf_effh_h
#define _chemistry_qc_scf_effh_h

#include <math/scmat/blkiter.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/scf/scf.h>

namespace sc {

class AccumEffectiveH: public BlockedSCElementOp2 {
  protected:
    SCF *scf_;
    double coef_[18];

    virtual void init() =0;
    
    // hindex is 0 for the closed and 1 for the open shell fock matrix
    // shelli and shellj are 0 for closed, 1 for open, and 2 for virtual
    int index(int hindex, int shelli, int shellj);

    // converts an occupation number to a shell number
    int shell(double);

    double& coef(int i, int j, int k) { return coef_[index(i,j,k)]; }

  public:
    AccumEffectiveH(SCF*);
    virtual ~AccumEffectiveH();

    virtual void process(SCMatrixBlockIter&,SCMatrixBlockIter&);
};

//  Guest & Saunders general form 
//        C        O         V
//    ----------
//    |        |
// C  |   fc   |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo |   fc   |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   |   fc   |
//    |        |        |        |
//    ----------------------------
class GSGeneralEffH: public AccumEffectiveH {
  protected:
    void init();
    
  public:
    GSGeneralEffH(SCF*);
    ~GSGeneralEffH();
};

//  Guest & Saunders' form for high spin
//        C        O         V
//    ----------
//    |        |
// C  | 2fc-fo |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo | 2fc-fo |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   | 2fc-fo |
//    |        |        |        |
//    ----------------------------
class GSHighSpinEffH: public AccumEffectiveH {
  protected:
    void init();

  public:
    GSHighSpinEffH(SCF*);
    ~GSHighSpinEffH();
};

//  test form
//        C        O         V
//    ----------
//    |        |
// C  |   fo   |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo |   fo   |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   |   fo   |
//    |        |        |        |
//    ----------------------------
class TestEffH: public AccumEffectiveH {
  protected:
    void init();

  public:
    TestEffH(SCF*);
    ~TestEffH();
};

//  form for converged wavefunction
//        C        O         V
//    ----------
//    |        |
// C  |   fc   |
//    |        |
//    -------------------
//    |        |        |
// O  | 2fc-fo |   fo   |
//    |        |        |
//    ----------------------------
//    |        |        |        |
// V  |   fc   |   fo   |   fo   |
//    |        |        |        |
//    ----------------------------
class PsiEffH: public AccumEffectiveH {
  protected:
    void init();

  public:
    PsiEffH(SCF*);
    ~PsiEffH();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
