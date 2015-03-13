//
// lgbuild.h --- definition of the local G matrix builder
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

#ifndef _chemistry_qc_scf_lgbuild_h
#define _chemistry_qc_scf_lgbuild_h

#undef SCF_CHECK_INTS
#undef SCF_CHECK_BOUNDS
#undef SCF_DONT_USE_BOUNDS

#include <mpqc_config.h>
#include <chemistry/qc/scf/gbuild.h>

namespace sc {

template<class T>
class LocalGBuild : public GBuild<T> {
  public:
    double tnint;    //!< total # of integrals computed in this thread
    double tinttime; //!< total user time spent computing integrals
    
  protected:
    MessageGrp *grp_;
    TwoBodyInt *tbi_;
    GaussianBasisSet *gbs_;
    PetiteList *rpl_;

    signed char * RESTRICT pmax;
    int threadno_;
    int nthread_;
    double accuracy_;
    
  public:
    LocalGBuild(T& t, const Ref<TwoBodyInt>& tbi, const Ref<PetiteList>& rpl,
                const Ref<GaussianBasisSet>& bs, const Ref<MessageGrp>& g,
                signed char *pm, double acc, int nt=1, int tn=0) :
      GBuild<T>(t),
      pmax(pm), threadno_(tn), nthread_(nt), accuracy_(acc)
    {
      grp_ = g.pointer();
      tbi_ = tbi.pointer();
      rpl_ = rpl.pointer();
      gbs_ = bs.pointer();
    }
    ~LocalGBuild() {}

    void run() {
      int tol = (int) (log(accuracy_)/log(2.0));
      int me=grp_->me();
      int nproc = grp_->n();
  
      // grab references for speed
      GaussianBasisSet& gbs = *gbs_;
      PetiteList& pl = *rpl_;
      TwoBodyInt& tbi = *tbi_;

      tbi.set_redundant(0);
      const double *intbuf = tbi.buffer();

      tnint=0;
      tinttime=0.0;
      sc_int_least64_t threadind=0;
      sc_int_least64_t ijklind=0;

      for (int i=0; i < gbs.nshell(); i++) {
        if (!pl.in_p1(i))
          continue;

        int fi=gbs.shell_to_function(i);
        int ni=gbs(i).nfunction();
        
        for (int j=0; j <= i; j++) {
          int oij = i_offset(i)+j;
          
          if (!pl.in_p2(oij))
            continue;

          int fj=gbs.shell_to_function(j);
          int nj=gbs(j).nfunction();
          int pmaxij = pmax[oij];

          for (int k=0; k <= i; k++, ijklind++) {
            if (ijklind%nproc != me)
              continue;

            threadind++;
            if (threadind % nthread_ != threadno_)
              continue;
            
            int fk=gbs.shell_to_function(k);
            int nk=gbs(k).nfunction();

            int pmaxijk=pmaxij, ptmp;
            if ((ptmp=pmax[i_offset(i)+k]-1) > pmaxijk) pmaxijk=ptmp;
            if ((ptmp=pmax[ij_offset(j,k)]-1) > pmaxijk) pmaxijk=ptmp;
        
            int okl = i_offset(k);
            for (int l=0; l <= (k==i?j:k); l++,okl++) {
              int pmaxijkl = pmaxijk;
              if ((ptmp=pmax[okl]) > pmaxijkl) pmaxijkl=ptmp;
              if ((ptmp=pmax[i_offset(i)+l]-1) > pmaxijkl) pmaxijkl=ptmp;
              if ((ptmp=pmax[ij_offset(j,l)]-1) > pmaxijkl) pmaxijkl=ptmp;

              int qijkl = pl.in_p4(oij,okl,i,j,k,l);
              if (!qijkl)
                continue;

#ifdef SCF_CHECK_BOUNDS
              double intbound = pow(2.0,double(tbi.log2_shell_bound(i,j,k,l)));
              double pbound   = pow(2.0,double(pmaxijkl));
              intbound *= qijkl;
              GBuild<T>::contribution.set_bound(intbound, pbound);
#else
#  ifndef SCF_DONT_USE_BOUNDS
              if (tbi.log2_shell_bound(i,j,k,l)+pmaxijkl < tol)
                continue;
#  endif
#endif

#ifdef MPQC_NEW_FEATURES
              const auto tstart = sc::RegionTimer::get_time_point();
#endif
              tbi.compute_shell(i,j,k,l);
#ifdef MPQC_NEW_FEATURES
              const auto tstop = RegionTimer::get_time_point();
              const std::chrono::duration<double> time_elapsed = tstop - tstart;
              tinttime += time_elapsed.count();
#endif

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

#ifdef SCF_CHECK_INTS
#ifdef HAVE_ISNAN
                      if (isnan(pki_int))
                        abort();
#endif
#endif
                      
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

            GBuild<T>::contribution.cont5(ij,kl,val);

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
            
            GBuild<T>::contribution.cont4(ij,kl,val);

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
            
            GBuild<T>::contribution.cont3(ij,kl,val);
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

          GBuild<T>::contribution.cont4(ij,kl,val);

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

          GBuild<T>::contribution.cont3(ij,kl,val);

        } else {
          /*
           * This integral contributes to J of G(ij)
           *
           * pkval = (ijkl)
           */
          ij = i_offset(ii)+jj;
          kl = i_offset(kk)+ll;
          val = (ij==kl) ? 0.5*pki_int : pki_int;

          GBuild<T>::contribution.cont1(ij,kl,val);

          /*
           * and to K1 of G(ik)
           *
           * pkval = -0.25 * (ijkl)
           */
          ij = ij_offset(ii,kk);
          kl = ij_offset(jj,ll);
          val = (ij==kl) ? 0.5*pki_int : pki_int;

          GBuild<T>::contribution.cont2(ij,kl,val);

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

            GBuild<T>::contribution.cont2(ij,kl,val);
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
          GBuild<T>::contribution.cont4(i_offset(ii)+jj,i_offset(kk)+ll,pki_int);

          /*
           * this integral also contributes to K1 and K2 of
           * G(il)
           *
           * pkval = -0.25 * ((ijkl)+(ikjl))
           *       = -0.5 * (ijkl)
           */
          GBuild<T>::contribution.cont3(ij_offset(ii,ll),ij_offset(kk,jj),pki_int);

        } else if (ii == kk || jj == ll) {
          /*
           * if i=k or j=l, then this integral contributes
           * to J and K2 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ilkj)
           *       = 0.75 * (ijkl)
           */
          GBuild<T>::contribution.cont4(i_offset(ii)+jj,i_offset(kk)+ll,pki_int);

          /*
           * this integral also contributes to K1 and K2 of
           * G(ik)
           *
           * pkval = -0.25 * ((ijkl)+(ilkj))
           *       = -0.5 * (ijkl)
           */
          GBuild<T>::contribution.cont3(ij_offset(ii,kk),ij_offset(jj,ll),pki_int);

        } else {
          /*
           * This integral contributes to J of G(ij)
           *
           * pkval = (ijkl)
           */
          GBuild<T>::contribution.cont1(i_offset(ii)+jj,i_offset(kk)+ll,pki_int);

          /*
           * and to K1 of G(ik)
           *
           * pkval = -0.25 * (ijkl)
           */
          GBuild<T>::contribution.cont2(ij_offset(ii,kk),ij_offset(jj,ll),pki_int);

          /*
           * and to K2 of G(il)
           *
           * pkval = -0.25 * (ijkl)
           */
          GBuild<T>::contribution.cont2(ij_offset(ii,ll),ij_offset(kk,jj),pki_int);
        }
      }
                    }
                  }
                }
              }

              tnint += (double) ni*nj*nk*nl;
            }
          }
        }
      }
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
