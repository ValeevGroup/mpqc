//
// tccont.h --- definition of the two-configuration fock contribution class
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

#ifndef _chemistry_qc_scf_tccont_h
#define _chemistry_qc_scf_tccont_h

#ifdef __GNUC__
#pragma interface
#endif

///////////////////////////////////////////////////////////////////////////

class LocalTCContribution {
  private:
    double * const gmata;
    double * const gmatb;
    double * const kmata;
    double * const kmatb;

    double * const pmata;
    double * const pmatb;
    double * const opmata;
    double * const opmatb;

  public:
    LocalTCContribution(double *ga, double *pa, double *gb, double *pb,
                        double *ka, double *opa, double *kb, double *opb) :
      gmata(ga), pmata(pa), gmatb(gb), pmatb(pb),
      kmata(ka), opmata(opa), kmatb(kb), opmatb(opb) {}
    ~LocalTCContribution() {}

    inline void cont1(int ij, int kl, double val) {
      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];
    }
    
    inline void cont2(int ij, int kl, double val) {
      val *= 0.25;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];

      kmata[ij] += val*opmata[kl];
      kmata[kl] += val*opmata[ij];

      kmatb[ij] += val*opmatb[kl];
      kmatb[kl] += val*opmatb[ij];
    }
    
    inline void cont3(int ij, int kl, double val) {
      val *= 0.5;
      gmata[ij] -= val*pmata[kl];
      gmata[kl] -= val*pmata[ij];

      gmatb[ij] -= val*pmatb[kl];
      gmatb[kl] -= val*pmatb[ij];

      kmata[ij] += val*opmata[kl];
      kmata[kl] += val*opmata[ij];

      kmatb[ij] += val*opmatb[kl];
      kmatb[kl] += val*opmatb[ij];
    }
    
    inline void cont4(int ij, int kl, double val) {
      gmata[ij] += 0.75*val*pmata[kl];
      gmata[kl] += 0.75*val*pmata[ij];

      gmatb[ij] += 0.75*val*pmatb[kl];
      gmatb[kl] += 0.75*val*pmatb[ij];

      kmata[ij] += 0.25*val*opmata[kl];
      kmata[kl] += 0.25*val*opmata[ij];

      kmatb[ij] += 0.25*val*opmatb[kl];
      kmatb[kl] += 0.25*val*opmatb[ij];
    }
    
    inline void cont5(int ij, int kl, double val) {
      val *= 0.5;
      gmata[ij] += val*pmata[kl];
      gmata[kl] += val*pmata[ij];

      gmatb[ij] += val*pmatb[kl];
      gmatb[kl] += val*pmatb[ij];

      kmata[ij] += val*opmata[kl];
      kmata[kl] += val*opmata[ij];

      kmatb[ij] += val*opmatb[kl];
      kmatb[kl] += val*opmatb[ij];
    }
};

class LocalTCGradContribution {
  private:
    double * const pmat;
    double * const pmata;
    double * const pmatb;
    double c1sq;
    double c2sq;
    double c1c2;

  public:
    LocalTCGradContribution(double *p, double *pa, double *pb,
                            double c1, double c2) :
      pmat(p), pmata(pa), pmatb(pb)
    {
      c1sq = c1*c1;
      c2sq = c2*c2;
      c1c2 = c1*c2;
    }
    ~LocalTCGradContribution() {}

    inline double cont1(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        c1sq*(pmata[ij]*pmat[kl] + pmat[ij]*pmata[kl]) +
        c2sq*(pmatb[ij]*pmat[kl] + pmat[ij]*pmatb[kl]) +
        0.5*c1sq*pmata[ij]*pmata[kl] +
        0.5*c2sq*pmatb[ij]*pmatb[kl];
    }

    inline double cont2(int ij, int kl) {
      return pmat[ij]*pmat[kl] +
        c1sq*(pmata[ij]*pmat[kl] + pmat[ij]*pmata[kl]) +
        c2sq*(pmatb[ij]*pmat[kl] + pmat[ij]*pmatb[kl]) -
        c1c2*(pmata[ij]*pmatb[kl] + pmatb[ij]*pmata[kl]);
    }
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
// End:
