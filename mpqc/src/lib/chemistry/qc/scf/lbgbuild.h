//
// lbgbuild.h --- definitino of the load-balanced local G matrix builder
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

#ifndef _chemistry_qc_scf_lbgbuild_h
#define _chemistry_qc_scf_lbgbuild_h

#ifdef __GNUC__
#pragma interface
#endif

#include <chemistry/qc/scf/gbuild.h>

namespace sc {

template<class T>
class LocalLBGBuild : public GBuild<T> {
  protected:
    Ref<MessageGrp> grp_;
    Ref<TwoBodyInt> tbi_;
    Ref<Integral> integral_;
    Ref<GaussianBasisSet> gbs_;
    signed char *pmax;
    
  public:
    LocalLBGBuild(T& t, const Ref<TwoBodyInt>& tbi, const Ref<Integral>& ints,
                const Ref<GaussianBasisSet>& bs, const Ref<MessageGrp>& g,
                signed char *pm) :
      GBuild<T>(t), grp_(g), tbi_(tbi), integral_(ints), gbs_(bs), pmax(pm) {}
    ~LocalLBGBuild() {}

    void build_gmat(double accuracy) {
      Timer tim("ao_gmat");
      tim.set_default("quartet");

      double tnint=0;
      int tol = (int) (log(accuracy)/log(2.0));
      int me=grp_->me();
      int nproc = grp_->n();
  
      Ref<PetiteList> rpl = integral_->petite_list();
  
      // grab references for speed
      GaussianBasisSet& gbs = *gbs_.pointer();
      PetiteList& pl = *rpl.pointer();
      TwoBodyInt& tbi = *tbi_.pointer();

      tbi.set_redundant(0);
      const double *intbuf = tbi.buffer();

      int inds[4];

      // node zero passes out indices
      if (me==0) {
        int i;
        for (i=0; i < gbs.nshell(); i++) {
          if (!pl.in_p1(i))
            continue;

          inds[0]=i;
          
          for (int j=0; j <= i; j++) {
            int oij = i_offset(i)+j;
          
            if (!pl.in_p2(oij))
              continue;
            
            inds[1]=j;

            tim.enter_default();
            int from;
            grp_->recvt(2323, &from, 1);
            grp_->sendt(from, 3232, inds, 4);
            tim.exit_default();
          }
        }

        // now clean up
        inds[0] = inds[1] = inds[2] = inds[3] = -1;
        
        for (i=1; i < nproc; i++) {
          int from;
          grp_->recvt(2323, &from, 1);
          grp_->sendt(from, 3232, inds, 4);
        }
      }
      
      // all other nodes do the work
      else {
        do {
          grp_->sendt(0, 2323, &me, 1);
          grp_->recvt(3232, inds, 4);

          int i=inds[0];
          int j=inds[1];
          
          if (i < 0)
            break;
          
          int fi=gbs.shell_to_function(i);
          int ni=gbs(i).nfunction();
        
          int oij = i_offset(i)+j;
          
          int fj=gbs.shell_to_function(j);
          int nj=gbs(j).nfunction();
          int pmaxij = pmax[oij];

          for (int k=0; k <= i; k++) {
          int fk=gbs.shell_to_function(k);
          int nk=gbs(k).nfunction();

          int pmaxijk=pmaxij, ptmp;
          if ((ptmp=pmax[i_offset(i)+k]-2) > pmaxijk) pmaxijk=ptmp;
          if ((ptmp=pmax[ij_offset(j,k)]-2) > pmaxijk) pmaxijk=ptmp;
        
          int okl = i_offset(k);
          for (int l=0; l <= (k==i?j:k); l++,okl++) {
            
            int pmaxijkl = pmaxijk;
            if ((ptmp=pmax[okl]) > pmaxijkl) pmaxijkl=ptmp;
            if ((ptmp=pmax[i_offset(i)+l]-2) > pmaxijkl) pmaxijkl=ptmp;
            if ((ptmp=pmax[ij_offset(j,l)]-2) > pmaxijkl) pmaxijkl=ptmp;
              

            if (tbi.log2_shell_bound(i,j,k,l)+pmaxijkl < tol)
              continue;

            int qijkl = pl.in_p4(oij,okl,i,j,k,l);
            if (!qijkl)
              continue;

            tim.enter_default();
            tbi.compute_shell(i,j,k,l);
            tim.exit_default();

            int e12 = (i==j);
            int e34 = (k==l);
            int e13e24 = (i==k) && (j==l);
            int e_any = e12||e34||e13e24;
    
            int fl=gbs.shell_to_function(l);
            int nl=gbs(l).nfunction();
     
            int ii,jj,kk,ll;
            int I,J,K,L;
            int index=0;

            for (I=0, ii=fi; I < ni; I++, ii++) {
              for (J=0, jj=fj; J <= (e12 ? I : nj-1); J++, jj++) {
                for (K=0, kk=fk; K <= (e13e24 ? I : nk-1); K++, kk++) {
           
                  int lend = (e34 ? ((e13e24)&&(K==I) ? J : K)
                              : ((e13e24)&&(K==I)) ? J : nl-1);
                  for (L=0, ll=fl; L <= lend; L++, ll++, index++) {
                    double pki_int = intbuf[index];

                    if ((pki_int>0?pki_int:-pki_int) < 1.0e-15)
                      continue;

                    if (qijkl > 1)
                      pki_int *= qijkl;

      if (e_any) {
        int ij,kl;
        double val;

        if (jj == kk) {
          /*
           * if i=j=k or j=k=l, then this integral contributes
           * to J, K1, and K2 of G(ij), so
           * pkval = (ijkl) - 0.25 * ((ikjl)-(ilkj))
           *       = 0.5 * (ijkl)
           */
          if (ii == jj || kk == ll) {
            ij = i_offset(ii)+jj;
            kl = i_offset(kk)+ll;
            val = (ij==kl) ? 0.5*pki_int : pki_int;

            contribution.cont5(ij,kl,val);

          } else {
            /*
             * if j=k, then this integral contributes
             * to J and K1 of G(ij)
             *
             * pkval = (ijkl) - 0.25 * (ikjl)
             *       = 0.75 * (ijkl)
             */
            ij = i_offset(ii)+jj;
            kl = i_offset(kk)+ll;
            val = (ij==kl) ? 0.5*pki_int : pki_int;
            
            contribution.cont4(ij,kl,val);

            /*
             * this integral also contributes to K1 and K2 of
             * G(il)
             *
             * pkval = -0.25 * ((ijkl)+(ikjl))
             *       = -0.5 * (ijkl)
             */
            ij = ij_offset(ii,ll);
            kl = ij_offset(kk,jj);
            val = (ij==kl) ? 0.5*pki_int : pki_int;
            
            contribution.cont3(ij,kl,val);
          }
        } else if (ii == kk || jj == ll) {
          /*
           * if i=k or j=l, then this integral contributes
           * to J and K2 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ilkj)
           *       = 0.75 * (ijkl)
           */
          ij = i_offset(ii)+jj;
          kl = i_offset(kk)+ll;
          val = (ij==kl) ? 0.5*pki_int : pki_int;

          contribution.cont4(ij,kl,val);

          /*
           * this integral also contributes to K1 and K2 of
           * G(ik)
           *
           * pkval = -0.25 * ((ijkl)+(ilkj))
           *       = -0.5 * (ijkl)
           */
          ij = ij_offset(ii,kk);
          kl = ij_offset(jj,ll);
          val = (ij==kl) ? 0.5*pki_int : pki_int;

          contribution.cont3(ij,kl,val);

        } else {
          /*
           * This integral contributes to J of G(ij)
           *
           * pkval = (ijkl)
           */
          ij = i_offset(ii)+jj;
          kl = i_offset(kk)+ll;
          val = (ij==kl) ? 0.5*pki_int : pki_int;

          contribution.cont1(ij,kl,val);

          /*
           * and to K1 of G(ik)
           *
           * pkval = -0.25 * (ijkl)
           */
          ij = ij_offset(ii,kk);
          kl = ij_offset(jj,ll);
          val = (ij==kl) ? 0.5*pki_int : pki_int;

          contribution.cont2(ij,kl,val);

          if ((ii != jj) && (kk != ll)) {
            /*
             * if i!=j and k!=l, then this integral also
             * contributes to K2 of G(il)
             *
             * pkval = -0.25 * (ijkl)
             *
             * note: if we get here, then ik can't equal jl,
             * so pkval wasn't multiplied by 0.5 above.
             */
            ij = ij_offset(ii,ll);
            kl = ij_offset(kk,jj);

            contribution.cont2(ij,kl,val);
          }
        }
      } else { // !e_any
        if (jj == kk) {
          /*
           * if j=k, then this integral contributes
           * to J and K1 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ikjl)
           *       = 0.75 * (ijkl)
           */
          contribution.cont4(i_offset(ii)+jj,i_offset(kk)+ll,pki_int);

          /*
           * this integral also contributes to K1 and K2 of
           * G(il)
           *
           * pkval = -0.25 * ((ijkl)+(ikjl))
           *       = -0.5 * (ijkl)
           */
          contribution.cont3(ij_offset(ii,ll),ij_offset(kk,jj),pki_int);

        } else if (ii == kk || jj == ll) {
          /*
           * if i=k or j=l, then this integral contributes
           * to J and K2 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ilkj)
           *       = 0.75 * (ijkl)
           */
          contribution.cont4(i_offset(ii)+jj,i_offset(kk)+ll,pki_int);

          /*
           * this integral also contributes to K1 and K2 of
           * G(ik)
           *
           * pkval = -0.25 * ((ijkl)+(ilkj))
           *       = -0.5 * (ijkl)
           */
          contribution.cont3(ij_offset(ii,kk),ij_offset(jj,ll),pki_int);

        } else {
          /*
           * This integral contributes to J of G(ij)
           *
           * pkval = (ijkl)
           */
          contribution.cont1(i_offset(ii)+jj,i_offset(kk)+ll,pki_int);

          /*
           * and to K1 of G(ik)
           *
           * pkval = -0.25 * (ijkl)
           */
          contribution.cont2(ij_offset(ii,kk),ij_offset(jj,ll),pki_int);

          /*
           * and to K2 of G(il)
           *
           * pkval = -0.25 * (ijkl)
           */
          contribution.cont2(ij_offset(ii,ll),ij_offset(kk,jj),pki_int);
        }
      }
                  }
                }
              }
            }

            tnint += (double) ni*nj*nk*nl;
          }
          }
        } while (inds[0] > -1);
      }

      grp_->sum(&tnint, 1, 0, 0);
      ExEnv::out0() << indent << scprintf("%20.0f integrals\n", tnint);

      tim.exit("ao_gmat");
    }
    
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
