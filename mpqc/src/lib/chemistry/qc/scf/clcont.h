//
// clcont.h --- definition of the closed-shell fock contribution class
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

#ifndef _chemistry_qc_scf_clcont_h
#define _chemistry_qc_scf_clcont_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/scf.h>

///////////////////////////////////////////////////////////////////////////

class LocalCLContribution {
  private:
    double * const gmat;
    double * const pmat;
    double bound;
  public:
    LocalCLContribution(double *g, double *p) : gmat(g), pmat(p), bound(0.0) {}
    ~LocalCLContribution() {}

    void set_bound(double b) { bound = b; }

    inline void cont1(int ij, int kl, double val) {
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
#if SCF_CHECK_BOUNDS
      if (fabs(val*pmat[kl]) > bound
          || fabs(val*pmat[ij]) > bound) {
        cerr << "val=" << val
             << ", pmat[" << kl << "]=" << pmat[kl]
             << ", pmat[" << ij << "]=" << pmat[ij]
             << ", bound=" << bound
             << endl;
        abort();
      }
#endif
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= -0.25;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
#if SCF_CHECK_BOUNDS
      if (fabs(val*pmat[kl]) > bound
          || fabs(val*pmat[ij]) > bound) {
        cerr << "val=" << val
             << ", pmat[" << kl << "]=" << pmat[kl]
             << ", pmat[" << ij << "]=" << pmat[ij]
             << ", bound=" << bound
             << endl;
        abort();
      }
#endif
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= -0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
#if SCF_CHECK_BOUNDS
      if (fabs(val*pmat[kl]) > bound
          || fabs(val*pmat[ij]) > bound) {
        cerr << "val=" << val
             << ", pmat[" << kl << "]=" << pmat[kl]
             << ", pmat[" << ij << "]=" << pmat[ij]
             << ", bound=" << bound
             << endl;
        abort();
      }
#endif
    }
    
    inline void cont4(int ij, int kl, double val) {
      val *= 0.75;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
#if SCF_CHECK_BOUNDS
      if (fabs(val*pmat[kl]) > bound
          || fabs(val*pmat[ij]) > bound) {
        cerr << "val=" << val
             << ", pmat[" << kl << "]=" << pmat[kl]
             << ", pmat[" << ij << "]=" << pmat[ij]
             << ", bound=" << bound
             << endl;
        abort();
      }
#endif
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmat[ij] += val*pmat[kl];
      gmat[kl] += val*pmat[ij];
#if SCF_CHECK_BOUNDS
      if (fabs(val*pmat[kl]) > bound
          || fabs(val*pmat[ij]) > bound) {
        cerr << "val=" << val
             << ", pmat[" << kl << "]=" << pmat[kl]
             << ", pmat[" << ij << "]=" << pmat[ij]
             << ", bound=" << bound
             << endl;
        abort();
      }
#endif
    }
};

class LocalCLGradContribution {
  private:
    double * const pmat;

  public:
    LocalCLGradContribution(double *p) : pmat(p) {}
    ~LocalCLGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl];
    }
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
