//
// ltbgrad.h --- definition of the local two-electron gradient builder
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

#ifndef _chemistry_qc_scf_ltbgrad_h
#define _chemistry_qc_scf_ltbgrad_h

#include <math.h>

#include <util/misc/regtime.h>
#include <math/scmat/offset.h>

#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/petite.h>

#include <chemistry/qc/scf/tbgrad.h>

namespace sc {
  
template<class T>
class LocalTBGrad : public TBGrad<T> {
  public:
    double *tbgrad;

  protected:
    MessageGrp *grp_;
    TwoBodyDerivInt *tbi_;
    GaussianBasisSet *gbs_;
    PetiteList *rpl_;
    Molecule *mol_;

    double pmax_;
    double accuracy_;

    int threadno_;
    int nthread_;

  public:
    LocalTBGrad(T& t, const Ref<TwoBodyDerivInt>& tbdi, const Ref<PetiteList>& pl,
                const Ref<GaussianBasisSet>& bs, const Ref<MessageGrp>& g,
                double *tbg, double pm, double a, int nt = 1, int tn = 0,
                double exchange_fraction = 1.0) :
      TBGrad<T>(t,exchange_fraction),
      tbgrad(tbg), pmax_(pm), accuracy_(a), threadno_(tn), nthread_(nt)
    {
      grp_ = g.pointer();
      gbs_ = bs.pointer();
      rpl_ = pl.pointer();
      tbi_ = tbdi.pointer();
      mol_ = gbs_->molecule().pointer();
    }

    ~LocalTBGrad() {}
    
    void run() {
      int me = grp_->me();
      int nproc = grp_->n();
      
      // grab ref for convenience
      GaussianBasisSet& gbs = *gbs_;
      Molecule& mol = *mol_;
      PetiteList& pl = *rpl_;
      TwoBodyDerivInt& tbi = *tbi_;
      
      // create vector to hold skeleton gradient
      double *tbint = new double[mol.natom()*3];
      memset(tbint, 0, sizeof(double)*mol.natom()*3);

      // for bounds checking
      int PPmax = (int) (log(6.0*pmax_*pmax_)/log(2.0));
      int threshold = (int) (log(accuracy_)/log(2.0));
  
      int kindex=0;
      int threadind=0;
      for (int i=0; i < gbs.nshell(); i++) {
        if (!pl.in_p1(i))
          continue;
    
        int ni=gbs(i).nfunction();
        int fi=gbs.shell_to_function(i);
    
        for (int j=0; j <= i; j++) {
          int ij=i_offset(i)+j;
          if (!pl.in_p2(ij))
            continue;
      
          if (tbi.log2_shell_bound(i,j,-1,-1)+PPmax < threshold)
            continue;
      
          int nj=gbs(j).nfunction();
          int fj=gbs.shell_to_function(j);
    
          for (int k=0; k <= i; k++,kindex++) {
            if (kindex%nproc != me)
              continue;
            
            threadind++;
            if (threadind % nthread_ != threadno_)
              continue;
            
            int nk=gbs(k).nfunction();
            int fk=gbs.shell_to_function(k);
    
            for (int l=0; l <= ((i==k)?j:k); l++) {
              if (tbi.log2_shell_bound(i,j,k,l)+PPmax < threshold)
                continue;
          
              int kl=i_offset(k)+l;
              int qijkl;
              if (!(qijkl=pl.in_p4(ij,kl,i,j,k,l)))
                continue;
          
              int nl=gbs(l).nfunction();
              int fl=gbs.shell_to_function(l);

              DerivCenters cent;
              tbi.compute_shell(i,j,k,l,cent);

              const double * buf = tbi.buffer();
          
              double cscl, escl;

              this->set_scale(cscl, escl, i, j, k, l);

              int indijkl=0;
              int nx=cent.n();
              //if (cent.has_omitted_center()) nx--;
              for (int x=0; x < nx; x++) {
                int ix=cent.atom(x);
                int io=cent.omitted_atom();
                for (int ixyz=0; ixyz < 3; ixyz++) {
                  double tx = tbint[ixyz+ix*3];
                  double to = tbint[ixyz+io*3];
                  
                  for (int ip=0, ii=fi; ip < ni; ip++, ii++) {
                    for (int jp=0, jj=fj; jp < nj; jp++, jj++) {
                      for (int kp=0, kk=fk; kp < nk; kp++, kk++) {
                        for (int lp=0, ll=fl; lp < nl; lp++, ll++, indijkl++) {
                          double contrib;
                          double qint = buf[indijkl]*qijkl;

                          contrib = cscl*qint*
                            TBGrad<T>::contribution.cont1(ij_offset(ii,jj),
                                               ij_offset(kk,ll));

                          tx += contrib;
                          to -= contrib;

                          contrib = escl*qint*
                            TBGrad<T>::contribution.cont2(ij_offset(ii,kk),
                                               ij_offset(jj,ll));

                          tx += contrib;
                          to -= contrib;

                          if (i!=j && k!=l) {
                            contrib = escl*qint*
                              TBGrad<T>::contribution.cont2(ij_offset(ii,ll),
                                                 ij_offset(jj,kk));

                            tx += contrib;
                            to -= contrib;
                          }
                        }
                      }
                    }
                  }

                  tbint[ixyz+ix*3] = tx;
                  tbint[ixyz+io*3] = to;
                }
              }
            }
          }
        }
      }
      
      CharacterTable ct = mol.point_group()->char_table();
      SymmetryOperation so;

      for (int alpha=0; alpha < mol.natom(); alpha++) {
        double tbx = tbint[alpha*3+0];
        double tby = tbint[alpha*3+1];
        double tbz = tbint[alpha*3+2];
        
        for (int g=1; g < ct.order(); g++) {
          so = ct.symm_operation(g);
          int ap = pl.atom_map(alpha,g);

          tbx += tbint[ap*3+0]*so(0,0) + tbint[ap*3+1]*so(1,0) +
                 tbint[ap*3+2]*so(2,0);
          tby += tbint[ap*3+0]*so(0,1) + tbint[ap*3+1]*so(1,1) +
                 tbint[ap*3+2]*so(2,1);
          tbz += tbint[ap*3+0]*so(0,2) + tbint[ap*3+1]*so(1,2) +
                 tbint[ap*3+2]*so(2,2);
        }
        double scl = 1.0/(double)ct.order();
        tbgrad[alpha*3+0] += tbx*scl;
        tbgrad[alpha*3+1] += tby*scl;
        tbgrad[alpha*3+2] += tbz*scl;
      }
    
      delete[] tbint;
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
