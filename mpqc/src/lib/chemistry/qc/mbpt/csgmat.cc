//
// csgmat.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
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

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/local.h>
#include <math/scmat/repl.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/scf/lgbuild.h>
#include <chemistry/qc/scf/clhftmpl.h>

using namespace sc;

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

#define INT_MAX1(n1) ((n1)-1)
#define INT_MAX2(e12,i,n2) ((e12)?(i):((n2)-1))
#define INT_MAX3(e13e24,i,n3) ((e13e24)?(i):((n3)-1))
#define INT_MAX4(e13e24,e34,i,j,k,n4) \
  ((e34)?(((e13e24)&&((k)==(i)))?(j):(k)) \
        :((e13e24)&&((k)==(i)))?(j):(n4)-1)

enum Access { Read, Write, Accum };
static RefSymmSCMatrix
get_local_data(const RefSymmSCMatrix& m, double*& p,
               const Ref<MessageGrp> &msg, Access access)
{
  RefSymmSCMatrix l = m;
  
  if (!dynamic_cast<LocalSymmSCMatrix*>(l.pointer())
      && !dynamic_cast<ReplSymmSCMatrix*>(l.pointer())) {
    Ref<SCMatrixKit> k = new ReplSCMatrixKit;
    l = k->symmmatrix(m.dim());
    l->convert(m);

    if (access == Accum)
      l->assign(0.0);
  } else if (msg->n() > 1 && access==Accum) {
    l = m.clone();
    l.assign(0.0);
  }

  if (dynamic_cast<ReplSymmSCMatrix*>(l.pointer()))
    p = dynamic_cast<ReplSymmSCMatrix*>(l.pointer())->get_data();
  else
    p = dynamic_cast<LocalSymmSCMatrix*>(l.pointer())->get_data();

  return l;
}

static signed char *
init_pmax(double *pmat_data, const Ref<GaussianBasisSet> &basis)
{
  double l2inv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);
  
  GaussianBasisSet& gbs = *basis.pointer();
  
  signed char * pmax = new signed char[ioff(gbs.nshell())];

  int ish, jsh, ij;
  for (ish=ij=0; ish < gbs.nshell(); ish++) {
    int istart = gbs.shell_to_function(ish);
    int iend = istart + gbs(ish).nfunction();
    
    for (jsh=0; jsh <= ish; jsh++,ij++) {
      int jstart = gbs.shell_to_function(jsh);
      int jend = jstart + gbs(jsh).nfunction();
      
      double maxp=0, tmp;

      for (int i=istart; i < iend; i++) {
        int ijoff = ioff(i) + jstart;
        for (int j=jstart; j < ((ish==jsh) ? i+1 : jend); j++,ijoff++)
          if ((tmp=::fabs(pmat_data[ijoff])) > maxp)
            maxp=tmp;
      }

      if (maxp <= tol)
        maxp=tol;

      long power = long(log(maxp)*l2inv);
      if (power < SCHAR_MIN) pmax[ij] = SCHAR_MIN;
      else if (power > SCHAR_MAX) pmax[ij] = SCHAR_MAX;
      else pmax[ij] = (signed char) power;
    }
  }

  return pmax;
}

/**************************************************************************
 *
 * calculate the closed shell G matrix
 * assume all matrices are held locally -- IMBN
 *
 * input:
 *   Gmat     = matrix containing old G matrix
 *   DPmat    = matrix containing density diff matrix
 *
 * on return:
 *   Gmat contains the new G matrix
 *
 * return 0 on success and -1 on failure
 */

int
MBPT2::make_cs_gmat_new(RefSymmSCMatrix& Gmat,
                        const RefSymmSCMatrix& DPmat)
{
  int i;
  int nthread = thr_->nthread();

  Timer tim("gmat");

  Ref<PetiteList> pl = integral()->petite_list(basis());

  Gmat.assign(0.0); 

  // scale the off-diagonal elements of DPmat by 2.0
  DPmat->scale(2.0);
  DPmat->scale_diagonal(0.5);

  // now try to figure out the matrix specialization we're dealing with
  // if we're using Local matrices, then there's just one subblock, or
  // see if we can convert G and P to local matrices

  if (debug_>1) {
    DPmat.print("DPmat before build");
  }

  // grab the data pointers from the G and P matrices
  double *gmat, *pmat;
  RefSymmSCMatrix gtmp = get_local_data(Gmat, gmat, msg_, Accum);
  RefSymmSCMatrix ptmp = get_local_data(DPmat, pmat, msg_, Read);

  signed char * pmax = init_pmax(pmat, basis());
  
  LocalGBuild<LocalCLHFContribution> **gblds =
    new LocalGBuild<LocalCLHFContribution>*[nthread];
  LocalCLHFContribution **conts = new LocalCLHFContribution*[nthread];

  double **gmats = new double*[nthread];
  gmats[0] = gmat;
    
  Ref<GaussianBasisSet> bs = basis();
  int ntri = ioff(bs->nbasis());

  for (i=0; i < nthread; i++) {
    if (i) {
      gmats[i] = new double[ntri];
      memset(gmats[i], 0, sizeof(double)*ntri);
      }
    conts[i] = new LocalCLHFContribution(gmats[i], pmat);
    gblds[i] = new LocalGBuild<LocalCLHFContribution>(*conts[i], tbints_[i],
         pl, bs, msg_, pmax, cphf_epsilon_/1000.0, nthread, i
      );
    
    thr_->add_thread(i, gblds[i]);
    }

  if (thr_->start_threads() < 0) {
    ExEnv::err0() << indent
         << "MBPT: csgmat: error starting threads" << std::endl;
    abort();
    }

  if (thr_->wait_threads() < 0) {
    ExEnv::err0() << indent
         << "MBPT: csgmat: error waiting for threads" << std::endl;
    abort();
    }
      
  double tnint=0;
  for (i=0; i < nthread; i++) {
    tnint += gblds[i]->tnint;

    if (i) {
      for (int j=0; j < ntri; j++)
        gmat[j] += gmats[i][j];
      delete[] gmats[i];
      }

    delete gblds[i];
    delete conts[i];
    }

  delete[] gmats;
  delete[] gblds;
  delete[] conts;
  delete[] pmax;
      
  msg_->sum(&tnint, 1, 0, 0);
  //ExEnv::out0() << indent << scprintf("%20.0f integrals\n", tnint);

  // if we're running on multiple processors, then sum the G matrix
  if (msg_->n() > 1)
      msg_->sum(gmat, ioff(basis()->nbasis()));

  // if we're running on multiple processors, or we don't have local
  // matrices, then accumulate gtmp back into G
  int local = (dynamic_cast<LocalSCMatrixKit*>(basis()->matrixkit().pointer()) ||
            dynamic_cast<ReplSCMatrixKit*>(basis()->matrixkit().pointer())) ? 1:0;
  if (!local || msg_->n() > 1)
    Gmat->convert_accumulate(gtmp);

  // now symmetrize the skeleton G matrix, placing the result in dd
  Gmat.scale(1.0/(double)pl->order());
  RefSymmSCMatrix Gmat_so(so_dimension(), basis_matrixkit());
  if (debug_>1) {
    Gmat.print("skeleton Gmat before symmetrize");
  }
  pl->symmetrize(Gmat,Gmat_so);
  if (debug_>1) {
    Gmat_so.print("Gmat in SO basis");
  }
  Gmat = pl->to_AO_basis(Gmat_so);
  if (debug_>1) {
    Gmat.print("Gmat in AO basis");
  }
  BlockedSymmSCMatrix *blocked_Gmat = dynamic_cast<BlockedSymmSCMatrix*>(Gmat.pointer());
  if (!blocked_Gmat || blocked_Gmat->nblocks() != 1) {
    ExEnv::outn() << "csgmat.cc: Gmat is wrong type" << std::endl;
    abort();
    }
  Gmat = blocked_Gmat->block(0);

  tim.exit("gmat");

  return 0;
}

int
MBPT2::make_cs_gmat(RefSymmSCMatrix& Gmat, double *DPmat)
{
  int errcod;

  Timer tim("gmat");

  errcod = make_g_d_nor(Gmat, DPmat, intbuf_);

  if (errcod != 0) {
    fprintf(stderr,"mbpt_gmat: trouble forming gmat 3\n");
    return -1;
    }

  tim.exit("gmat");

  return 0;
}

/************************************************************************
 *
 * Form the vector maxp; each element of maxp is the 2-based log of the
 * largest element (absolute value) in a block of the density matrix
 * (DPmat). The density matrix is of dimension nbasis x nbasis
 *
 ************************************************************************/

void
MBPT2::form_max_dens(double *DPmat, signed char *maxp)
{

  int i, j, k, l, ij;
  int isize, jsize, ioffset, joffset;
  double linv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);
  double ftmp, tmp;
  double *dpmat_ptr;

  for (i=0; i<basis()->nshell(); i++) {
    isize = basis()->shell(i).nfunction();
    ioffset = basis()->shell_to_function(i);
    for (j=0; j<=i; j++) {
      jsize = basis()->shell(j).nfunction();
      joffset = basis()->shell_to_function(j);
      tmp = 0.0;
      for (k=0; k<isize; k++) {
        dpmat_ptr = &DPmat[nbasis*(ioffset+k) + joffset];
        for (l=0; l<jsize; l++) {
          ftmp = fabs(*dpmat_ptr++);
          if (ftmp > tmp) tmp = ftmp;
          }
        }
      tmp = (tmp > tol) ? tmp : tol;
      ij = i*(i+1)/2 +j;
      maxp[ij] = (signed char) (log(tmp)*linv);
             /* log(tmp)/linv equals the 2-based log of tmp */
      }
    }

}

int
MBPT2::init_cs_gmat()
{
  tbint_ = integral()->electron_repulsion();
  tbint_->set_redundant(0);
  intbuf_ = tbint_->buffer();
  return 1;
}

void
MBPT2::done_cs_gmat()
{
  tbint_ = 0;
  intbuf_ = 0;
}

int
MBPT2::make_g_d_nor(RefSymmSCMatrix& Gmat,
                    double *DPmat, const double *mgdbuff)
{
  int tmax,imax,cpmax,pmaxijk=0;
  int pmaxik,pmaxjk,pmaxij=0;
  int i,j,k,l;
  int ij,kl;
  int n1,n2,n3,n4;
  int e12,e34,e13e24,e_any;
  int bf1,bf2,bf3,bf4;
  int i1,j1,k1,l1;
  int i2,j2,k2,l2;
  int ii,jj,kk,ll;
  int ij1;
  int lij,lkl;
  int index;
  int int_index,kindex;
  int nproc=msg_->n();
  int me=msg_->me();
  int s1,s2,s3,s4;
  int nbatri = (nbasis*(nbasis+1))/2;

  double tol = desired_gradient_accuracy() / 1000.0;
  if (min_orthog_res() < 1.0) { tol *= min_orthog_res(); }
  int inttol = (int) (log(tol)/log(2.0));

  double tnint=0.0;
  double pki_int,value;
  double *gtmp=0, *ptmp=0;
  double *dpmat_ptr;

  char *shnfunc=0;
  signed char *maxp=0;


  // Scale DPmat; this is necessary when using the gmat formation
  // program from scf (modified slightly), since this program assumes
  // that the off-diagonal elements have been scaled by a factor of 2.0
  dpmat_ptr = DPmat;
  for (i=0; i<nbasis; i++) {
    for (j=0; j<nbasis; j++) {
      if (i != j) *dpmat_ptr++ *= 2.0;
      else dpmat_ptr++;
      }
    }

  // Allocate and assign maxp
  if (eliminate_in_gmat_) {
    int nshellt = basis()->nshell()*(basis()->nshell()+1)/2;
    maxp = (signed char*) malloc(sizeof(signed char)*nshellt);
    if (!(maxp)) {
      fprintf(stderr,"mkgdlb: could not malloc maxp\n");
      return -1;
      }
    form_max_dens(DPmat, maxp);
    }

  // Allocate and assign ptmp (contains lower triangle of DPmat
  ptmp = (double*) malloc(sizeof(double)*nbatri);
  if (!(ptmp)) {
    fprintf(stderr,"mkgdlb: could not malloc ptmp\n");
    return -1;
    }
  for (i=0; i<nbasis; i++) {
    dpmat_ptr = &DPmat[i*nbasis];
    for (j=0; j<=i; j++) {
      ptmp[i*(i+1)/2 + j] = *dpmat_ptr++;
      }
    }


  // "Unscale" DPmat to get the original DPmat
  dpmat_ptr = DPmat;
  for (i=0; i<nbasis; i++) {
    for (j=0; j<nbasis; j++) {
      if (i != j)  *dpmat_ptr++ *= 0.50;
      else dpmat_ptr++;
      }
    }

  // Allocate and initialize gtmp
  gtmp = (double *) malloc(sizeof(double)*nbatri);
  for (i=0; i<nbatri; i++) gtmp[i] = 0.0;
 
  // Allocate and assign shnfunc
  shnfunc = (char *) malloc(basis()->nshell());
  if (!shnfunc) {
    fprintf(stderr,"make_g_d_lb: could not malloc shnfunc\n");
    return -1;
    }
  for (i=0; i < basis()->nshell(); i++) shnfunc[i]=basis()->shell(i).nfunction();


  /********************************************************
   * Start the actual formation of the G matrix:          *
   * Loop over all shells, calculate a bunch of integrals *
   * from each shell quartet, and stick those integrals   *
   * where they belong                                    *
   ********************************************************/


  kindex=int_index=0;
  for (i=0; i<basis()->nshell(); i++) {

    for (j=0; j<=i; j++) {
      ij = ioff(i)+j;
      if(eliminate_in_gmat_) pmaxij=maxp[ij];

      for (k=0; k<=i; k++,kindex++) {
        if(kindex%nproc!=me) {
          continue;
          }

        kl=ioff(k);
        if(eliminate_in_gmat_) {
          pmaxijk=pmaxij;
          if((pmaxik=maxp[(ioff(i)+k)]-2)>pmaxijk) pmaxijk=pmaxik;
          if((pmaxjk=maxp[IOFF(j,k)]-2)>pmaxijk) pmaxijk=pmaxjk;
          }

        for (l=0; l<=(k==i?j:k); l++) {

          imax = (int) tbint_->log2_shell_bound(i,j,k,l);

          if(eliminate_in_gmat_) {
            cpmax = (maxp[kl]>pmaxijk) ? maxp[kl] : pmaxijk;
            if((tmax=maxp[(ioff(i)+l)]-2)>cpmax) cpmax=tmax;
            if((tmax=maxp[IOFF(j,l)]-2)>cpmax) cpmax=tmax;

            if(cpmax+imax < inttol) {
              kl++;
              continue;
              }
            }

            s1 = i; s2 = j; s3 = k; s4 = l;

            tbint_->compute_shell(s1,s2,s3,s4);

            n1 = shnfunc[s1];
            n2 = shnfunc[s2];
            n3 = shnfunc[s3];
            n4 = shnfunc[s4];

           // Shell equivalence information
            e12    = (s2==s1);
            e13e24 = (s3==s1) && (s4==s2);
            e34    = (s4==s3);

            index = 0;

            e_any = (e12||e13e24||e34);
            if(e_any) {
              for (bf1=0; bf1<=INT_MAX1(n1) ; bf1++) {
                i2 = basis()->shell_to_function(s1) + bf1;

                for (bf2=0; bf2<=INT_MAX2(e12,bf1,n2) ; bf2++) {
                  j2 = basis()->shell_to_function(s2) + bf2;
                  if(i2>=j2) { i1=i2; j1=j2; }
                  else { i1=j2; j1=i2; }
                  ij1=ioff(i1)+j1;

                  for (bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3) ; bf3++) {
                    k2 = basis()->shell_to_function(s3) + bf3;

                    for (bf4=0;bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4);bf4++){
                      if (fabs(mgdbuff[index])>1.0e-10) {
                        l2 = basis()->shell_to_function(s4) + bf4;

                        if(k2>=l2) { k1=k2; l1=l2; }
                        else { k1=l2; l1=k2; }

                        if(ij1 >= ioff(k1)+l1) {
                          ii = i1; jj = j1; kk = k1; ll = l1;
                          }
                        else {
                          ii = k1; jj = l1; kk = i1; ll = j1;
                          }

                        pki_int = mgdbuff[index];

                        if (jj == kk) {
                          if (ii == jj || kk == ll) {
                            lij=ioff(ii)+jj;
                            lkl=ioff(kk)+ll;
                            value=(lij==lkl)? 0.25*pki_int: 0.5*pki_int;
                            gtmp[lij] += ptmp[lkl]*value;
                            gtmp[lkl] += ptmp[lij]*value;
                            }
                          else {
                            lij=ioff(ii)+jj;
                            lkl=ioff(kk)+ll;
                            value=(lij==lkl)? 0.375*pki_int: 0.75*pki_int;
                            gtmp[lij] += ptmp[lkl]*value;
                            gtmp[lkl] += ptmp[lij]*value;

                            lij=ioff(ii)+ll;
                            lkl=IOFF(kk,jj);
                            value=(lij==lkl)? 0.25*pki_int: 0.5*pki_int;
                            gtmp[lij] -= ptmp[lkl]*value;
                            gtmp[lkl] -= ptmp[lij]*value;
                            }
                          }
                        else if (ii == kk || jj == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=(lij==lkl)? 0.375*pki_int: 0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=(lij==lkl)? 0.25*pki_int : 0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        else {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=(lij==lkl)? 0.5*pki_int : pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=(lij==lkl)? 0.125*pki_int: 0.25*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;

                          if((ii != jj) && (kk != ll)) {
                            lij=ioff(ii)+ll;
                            lkl=IOFF(kk,jj);
                            value=(lij==lkl)? 0.125*pki_int: 0.25*pki_int;
                            gtmp[lij] -= ptmp[lkl]*value;
                            gtmp[lkl] -= ptmp[lij]*value;
                            }
                          }
                        }
                      index++;
                      }
                    }
                  }
                }
              }
            else {
              for (bf1=0; bf1<n1 ; bf1++) {
                i2 = basis()->shell_to_function(s1) + bf1;

                for (bf2=0; bf2<n2 ; bf2++) {
                  j2 = basis()->shell_to_function(s2) + bf2;
                  if(i2>=j2) { i1=i2; j1=j2; }
                  else { i1=j2; j1=i2; }
                  ij1=ioff(i1)+j1;

                  for (bf3=0; bf3<n3 ; bf3++) {
                    k2 = basis()->shell_to_function(s3) + bf3;

                    for (bf4=0; bf4<n4; bf4++) {
                      if (fabs(mgdbuff[index])>1.0e-10) {
                        l2 = basis()->shell_to_function(s4) + bf4;

                        if(k2>=l2) { k1=k2; l1=l2; }
                        else { k1=l2; l1=k2; }

                        if(ij1 >= ioff(k1)+l1) {
                          ii = i1; jj = j1; kk = k1; ll = l1;
                          }
                        else {
                          ii = k1; jj = l1; kk = i1; ll = j1;
                          }

                        pki_int = mgdbuff[index];

                        if (jj == kk) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          value=0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        else if (ii == kk || jj == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        else {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
  
                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value*=0.25;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;

                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          }
                        }
                      index++;
                      }
                    }
                  }
                }
              }
            tnint += (double) (n1*n2*n3*n4);
          kl++;
          int_index++;
          } // exit l loop
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop

  // Sum up contributions to gtmp
  msg_->sum(gtmp,nbatri,ptmp);


  // Put gtmp back into Gmat
  for (i=0; i<nbasis; i++) {
    for (j=0; j<=i; j++) {
      ij = i*(i+1)/2 + j;
      Gmat->set_element(i,j,gtmp[ij]);
//    Gmat->set_element(j,i,gtmp[ij]);  don't do this - only lower triangle
      }
    }


  // Free up memory
  if (gtmp) free(gtmp);
  if (maxp) free(maxp);
  if (ptmp) free(ptmp);
  if (shnfunc) free(shnfunc);

  return 0;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
