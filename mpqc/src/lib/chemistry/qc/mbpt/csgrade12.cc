//
// csgrade12.cc
// based on: csgrad.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/csgrade12.h>

#define PRINT1Q 0

/////////////////////////////////////////////////////////////////
// Function iquicksort performs a quick sort (larger -> smaller) 
// of the integer data in item by the integer indices in index;
// data in item remain unchanged
/////////////////////////////////////////////////////////////////
static void
iqs(int *item,int *index,int left,int right)
{
  register int i,j;
  int x,y;
 
  i=left; j=right;
  x=item[index[(left+right)/2]];
 
  do {
    while(item[index[i]]>x && i<right) i++;
    while(x>item[index[j]] && j>left) j--;
 
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);
       
  if (left<j) iqs(item,index,left,j);
  if (i<right) iqs(item,index,i,right);
}

static void
iquicksort(int *item,int *index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  iqs(item,index,0,n-1);
  }

CSGradErep12Qtr::CSGradErep12Qtr(int mythread_a, int nthread_a,
                                 int me_a, int nproc_a,
                                 const RefMemoryGrp &mem_a,
                                 const RefMessageGrp &msg_a,
                                 const RefThreadLock &lock_a,
                                 const RefGaussianBasisSet &basis_a,
                                 const RefTwoBodyInt &tbint_a,
                                 int ni_a, int nocc_a,
                                 double **scf_vector_a,
                                 double tol_a, int debug_a,
                                 int dynamic_a)
{
  msg = msg_a;
  mythread = mythread_a;
  nthread = nthread_a;
  lock = lock_a;
  basis = basis_a;
  tbint = tbint_a;
  ni = ni_a;
  nocc = nocc_a;
  me = me_a;
  nproc = nproc_a;
  tol = tol_a;
  mem = mem_a;
  scf_vector = scf_vector_a;
  debug = debug_a;
  dynamic_ = dynamic_a;

  aoint_computed = 0;
  timer = new RegionTimer();
}

CSGradErep12Qtr::~CSGradErep12Qtr()
{
}

void
CSGradErep12Qtr::run()
{
  if (dynamic_ && nproc > 1) {
    run_dynamic();
    }
  else {
    run_static();
    }
}

void
CSGradErep12Qtr::run_dynamic()
{
  if (me == 0 && mythread == 0) {
    run_task_manager();
    }
  else if (me > 0) {
    run_task_runner();
    }
}

void
CSGradErep12Qtr::run_task_manager()
{
  // intialize work arrays
  int S,R,index;
  int nshell = basis->nshell();
  int ntri = (nshell*(nshell+1))/2;
  int *cost = new int[ntri];
  int *Svec = new int[ntri];
  int *Rvec = new int[ntri];
  int *Ivec = new int[ntri];
  index = 0;
  for (S=0; S<nshell; S++) {
    for (R=0; R<=S; R++) {
      cost[index] = basis->shell(S).nfunction()*basis->shell(R).nfunction();
      Svec[index] = S;
      Rvec[index] = R;
      Ivec[index] = index;
      index++;
      }
    }

  // sort work
  iquicksort(cost, Ivec, ntri);
  if (debug) {
    cout << "costs of shell pairs" << endl;
    for (index=0; index<ntri; index++) {
      cout << scprintf(" (%d %d):%d",Svec[Ivec[index]],Rvec[Ivec[index]],
                       cost[Ivec[index]])
           << endl;
      }
    }

  // process requests
  int nreq = ntri + nthread*(nproc - 1);
  int iwork = 0;
  int print_index = 0;
  int print_interval = nreq/100;
  if (print_interval==0) print_interval = 1;
  int nreq_left = nreq;
  while (nreq_left) {
    int node;
    msg->recvt(18101,&node,1);
    int SR[2];
    if (iwork < ntri) {
      SR[0] = Svec[Ivec[iwork]];
      SR[1] = Rvec[Ivec[iwork]];
      iwork++;
      }
    else {
      SR[0] = -1;
      SR[1] = -1;
      }
    if (debug && print_index++%print_interval == 0) {
      cout << scprintf("sending %3d %3d to %3d, %3d%% complete",
                       SR[0],SR[1],node,(print_index*100)/nreq)
           << endl;
      }
    msg->sendt(node,18102,SR,2);
    nreq_left--;
    }

  if (debug) {
      cout << "all requests processed" << endl;
    }

  delete[] cost;
  delete[] Svec;
  delete[] Rvec;
  delete[] Ivec;
}

int
CSGradErep12Qtr::get_task(int &S, int &R)
{
  int SR[2];

  lock->lock();
  msg->sendt(0,18101,&me,1);
  msg->recvt(18102,SR,2);
  lock->unlock();

  S = SR[0];
  R = SR[1];
  if (S == -1) return 0;
  return 1;
}

void
CSGradErep12Qtr::run_task_runner()
{
  int P,Q,R,S;
  int p,q,r,s;
  int np,nq,nr,ns;
  int bf1,bf2,bf3,bf4;
  int p_offset,q_offset,r_offset,s_offset;
  int index = 0;
  int thindex = 0;
  int offset;
  int nfuncmax = basis->max_nfunction_in_shell();
  int nshell = basis->nshell();
  int nbasis = basis->nbasis();
  double dtol = pow(2.0,tol);
  double *iqjs_ptr;
  double *iqrs_ptr, *iprs_ptr;
  double *c_pi, *c_qi;
  double tmpval;
  int i,j;
  double *iqjs_contrib;  // local contributions to integral_iqjs
  double *iqjr_contrib;  // local contributions to integral_iqjr
  const int catchup_mask = 3;
  int catchup_ctr = 0;

  const double *intbuf = tbint->buffer();

  iqjs_contrib  = new double[nbasis*nfuncmax];
  iqjr_contrib  = new double[nbasis*nfuncmax];

  double *integral_iqrs; // quarter transformed two-el integrals
  lock->lock();
  integral_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
  lock->unlock();

  int work_per_thread = ((nshell*(nshell+1))/2)/(nproc*nthread);
  int print_interval = work_per_thread/100;
  int time_interval = work_per_thread/10;
  int print_index = 0;
  if (print_interval == 0) print_interval = 1;
  if (time_interval == 0) time_interval = 1;
  if (work_per_thread == 0) work_per_thread = 1;

  if (debug) {
      cout << scprintf("%d:%d: starting get_task loop",me,mythread) << endl;
    }

  S = 0;
  R = 0;
  while (get_task(S,R)) {
    ns = basis->shell(S).nfunction();
    s_offset = basis->shell_to_function(S);

    nr = basis->shell(R).nfunction();
    r_offset = basis->shell_to_function(R);

    if (debug > 1 && (print_index++)%print_interval == 0) {
      lock->lock();
      cout << scprintf("%d:%d: (PQ|%d %d) %d%%",
                       me,mythread,R,S,(100*print_index)/work_per_thread)
           << endl;
      lock->unlock();
      }
    if (debug > 1 && (print_index)%time_interval == 0) {
      lock->lock();
      cout << scprintf("timer for %d:%d:",me,mythread) << endl;
      timer->print();
      lock->unlock();
      }

    bzerofast(integral_iqrs, ni*nbasis*nfuncmax*nfuncmax);

    for (Q=0; Q<nshell; Q++) {
      nq = basis->shell(Q).nfunction();
      q_offset = basis->shell_to_function(Q);
      for (P=0; P<=Q; P++) {
        np = basis->shell(P).nfunction();
        p_offset = basis->shell_to_function(P);

        if (tbint->log2_shell_bound(P,Q,R,S) < tol) {
          continue;  // skip ereps less than tol
          }

        aoint_computed++;

        timer->enter("erep");
        tbint->compute_shell(P,Q,R,S);
        timer->exit("erep");

        lock->lock(); mem->catchup(); lock->unlock();

        timer->enter("1. q.t.");
        // Begin first quarter transformation;
        // generate (iq|rs) for i active

        offset = nr*ns*nbasis;
        const double *pqrs_ptr = intbuf;
        for (bf1 = 0; bf1 < np; bf1++) {
          p = p_offset + bf1;
          for (bf2 = 0; bf2 < nq; bf2++) {
            q = q_offset + bf2;

            if (q < p) {
              pqrs_ptr = &intbuf[ns*nr*(bf2+1 + nq*bf1)];
              continue; // skip to next q value
              }

            for (bf3 = 0; bf3 < nr; bf3++) {
              r = r_offset + bf3;

              for (bf4 = 0; bf4 < ns; bf4++) {
                s = s_offset + bf4;

                if (s < r) {
                  pqrs_ptr++;
                  continue; // skip to next bf4 value
                  }

                if (fabs(*pqrs_ptr) > dtol) {
                  iprs_ptr = &integral_iqrs[bf4 + ns*(p + nbasis*bf3)];
                  iqrs_ptr = &integral_iqrs[bf4 + ns*(q + nbasis*bf3)];
                  c_qi = &scf_vector[q][i_offset];
                  c_pi = &scf_vector[p][i_offset];
                  tmpval = *pqrs_ptr;
                  for (i=0; i<ni; i++) {
                    *iprs_ptr += *c_qi++*tmpval;
                    iprs_ptr += offset;
                    if (p != q) {
                      *iqrs_ptr += *c_pi++*tmpval;
                      iqrs_ptr += offset;
                      }
                    } // exit i loop
                  }   // endif

                pqrs_ptr++;
                } // exit bf4 loop
              }   // exit bf3 loop
            }     // exit bf2 loop
          }       // exit bf1 loop
        // end of first quarter transformation
        timer->exit("1. q.t.");

        }           // exit P loop
      }             // exit Q loop

#if PRINT1Q
      {
      lock->lock();
      double *tmp = integral_iqrs;
      for (int i = 0; i<ni; i++) {
        for (int r = 0; r<nr; r++) {
          for (int q = 0; q<nbasis; q++) {
            for (int s = 0; s<ns; s++) {
              printf("1Q: (%d %d|%d %d) = %12.8f\n",
                     i,q,r+r_offset,s+s_offset,*tmp);
              tmp++;
              }
            }
          }
        }
      lock->unlock();
      }
#endif

    timer->enter("2. q.t.");
    // Begin second quarter transformation;
    // generate (iq|jr) for i active and j active or frozen
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {

        bzerofast(iqjs_contrib, nbasis*nfuncmax);
        bzerofast(iqjr_contrib, nbasis*nfuncmax);

        for (bf1=0; bf1<ns; bf1++) {
          s = s_offset + bf1;
          double *c_sj = &scf_vector[s][j];
          double *iqjr_ptr = iqjr_contrib;
          for (bf2=0; bf2<nr; bf2++) {
            r = r_offset + bf2;
            if (r > s) {
              break; // skip to next bf1 value
              }
            double c_rj = scf_vector[r][j];
            iqjs_ptr = &iqjs_contrib[bf1*nbasis];
            iqrs_ptr = &integral_iqrs[bf1 + ns*nbasis*(bf2 + nr*i)];
            for (q=0; q<nbasis; q++) {
              *iqjs_ptr++ += c_rj * *iqrs_ptr;
              if (r != s) *iqjr_ptr += *c_sj * *iqrs_ptr;
              iqjr_ptr++;
              iqrs_ptr += ns;
              } // exit q loop
            // every so often process outstanding messages
            if ((catchup_ctr++ & catchup_mask) == 0) { lock->lock(); mem->catchup(); lock->unlock(); }
            }   // exit bf2 loop
          }     // exit bf1 loop

        // We now have contributions to iqjs and iqjr for one pair i,j,
        // all q, r in R and s in S; send iqjs and iqjr to the node
        // (ij_proc) which is going to have this ij pair
        int ij_proc =  (i*nocc + j)%nproc;
        int ij_index = (i*nocc + j)/nproc;

        // Sum the iqjs_contrib to the appropriate place
        int ij_offset = nbasis*(s_offset + nbasis*ij_index);
        lock->lock();
        mem->sum_reduction_on_node(iqjs_contrib,
                                   ij_offset, ns*nbasis, ij_proc);

        ij_offset = nbasis*(r_offset + nbasis*ij_index);
        mem->sum_reduction_on_node(iqjr_contrib,
                                   ij_offset, nr*nbasis, ij_proc);
        lock->unlock();

        }     // exit j loop
      }       // exit i loop
    // end of second quarter transformation
    timer->exit("2. q.t.");

    }         // exit while get_task

  if (debug) {
      cout << scprintf("%d:%d: done with get_task loop",me,mythread) << endl;
    }

  lock->lock();
  delete[] integral_iqrs;
  delete[] iqjs_contrib;
  delete[] iqjr_contrib;
  lock->unlock();
}

void
CSGradErep12Qtr::run_static()
{
  int P,Q,R,S;
  int p,q,r,s;
  int np,nq,nr,ns;
  int bf1,bf2,bf3,bf4;
  int p_offset,q_offset,r_offset,s_offset;
  int index = 0;
  int thindex = 0;
  int offset;
  int nfuncmax = basis->max_nfunction_in_shell();
  int nshell = basis->nshell();
  int nbasis = basis->nbasis();
  double dtol = pow(2.0,tol);
  double *iqjs_ptr;
  double *iqrs_ptr, *iprs_ptr;
  double *c_pi, *c_qi;
  double tmpval;
  int i,j;
  double *iqjs_contrib;  // local contributions to integral_iqjs
  double *iqjr_contrib;  // local contributions to integral_iqjr
  const int catchup_mask = 3;
  int catchup_ctr = 0;

  const double *intbuf = tbint->buffer();

  iqjs_contrib  = new double[nbasis*nfuncmax];
  iqjr_contrib  = new double[nbasis*nfuncmax];

  double *integral_iqrs; // quarter transformed two-el integrals
  lock->lock();
  integral_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
  lock->unlock();

  int work_per_thread = ((nshell*(nshell+1))/2)/(nproc*nthread);
  int print_interval = work_per_thread/100;
  int time_interval = work_per_thread/10;
  int print_index = 0;
  if (print_interval == 0) print_interval = 1;
  if (time_interval == 0) time_interval = 1;
  if (work_per_thread == 0) work_per_thread = 1;

  for (S=0; S<nshell; S++) {
    ns = basis->shell(S).nfunction();
    s_offset = basis->shell_to_function(S);

    for (R=0; R<=S; R++) {
      nr = basis->shell(R).nfunction();
      r_offset = basis->shell_to_function(R);

      if (index++%nproc == me && thindex++%nthread == mythread) {

        if (debug && (print_index++)%print_interval == 0) {
          lock->lock();
          cout << scprintf("%d:%d: (PQ|%d %d) %d%%",
                           me,mythread,R,S,(100*print_index)/work_per_thread)
               << endl;
          lock->unlock();
          }
        if (debug && (print_index)%time_interval == 0) {
          lock->lock();
          cout << scprintf("timer for %d:%d:",me,mythread) << endl;
          timer->print();
          lock->unlock();
          }

        bzerofast(integral_iqrs, ni*nbasis*nfuncmax*nfuncmax);

        for (Q=0; Q<nshell; Q++) {
          nq = basis->shell(Q).nfunction();
          q_offset = basis->shell_to_function(Q);
          for (P=0; P<=Q; P++) {
            np = basis->shell(P).nfunction();
            p_offset = basis->shell_to_function(P);

            if (tbint->log2_shell_bound(P,Q,R,S) < tol) {
              continue;  // skip ereps less than tol
              }

            aoint_computed++;

            timer->enter("erep");
            tbint->compute_shell(P,Q,R,S);
            timer->exit("erep");

            lock->lock(); mem->catchup(); lock->unlock();

            timer->enter("1. q.t.");
            // Begin first quarter transformation;
            // generate (iq|rs) for i active

            offset = nr*ns*nbasis;
            const double *pqrs_ptr = intbuf;
            for (bf1 = 0; bf1 < np; bf1++) {
              p = p_offset + bf1;
              for (bf2 = 0; bf2 < nq; bf2++) {
                q = q_offset + bf2;

                if (q < p) {
                  pqrs_ptr = &intbuf[ns*nr*(bf2+1 + nq*bf1)];
                  continue; // skip to next q value
                  }

                for (bf3 = 0; bf3 < nr; bf3++) {
                  r = r_offset + bf3;

                  for (bf4 = 0; bf4 < ns; bf4++) {
                    s = s_offset + bf4;

                    if (s < r) {
                      pqrs_ptr++;
                      continue; // skip to next bf4 value
                      }

                    if (fabs(*pqrs_ptr) > dtol) {
                      iprs_ptr = &integral_iqrs[bf4 + ns*(p + nbasis*bf3)];
                      iqrs_ptr = &integral_iqrs[bf4 + ns*(q + nbasis*bf3)];
                      c_qi = &scf_vector[q][i_offset];
                      c_pi = &scf_vector[p][i_offset];
                      tmpval = *pqrs_ptr;
                      for (i=0; i<ni; i++) {
                        *iprs_ptr += *c_qi++*tmpval;
                        iprs_ptr += offset;
                        if (p != q) {
                          *iqrs_ptr += *c_pi++*tmpval;
                          iqrs_ptr += offset;
                          }
                        } // exit i loop
                      }   // endif

                    pqrs_ptr++;
                    } // exit bf4 loop
                  }   // exit bf3 loop
                }     // exit bf2 loop
              }       // exit bf1 loop
            // end of first quarter transformation
            timer->exit("1. q.t.");

            }           // exit P loop
          }             // exit Q loop

#if PRINT1Q
          {
          lock->lock();
          double *tmp = integral_iqrs;
          for (int i = 0; i<ni; i++) {
            for (int r = 0; r<nr; r++) {
              for (int q = 0; q<nbasis; q++) {
                for (int s = 0; s<ns; s++) {
                  printf("1Q: (%d %d|%d %d) = %12.8f\n",
                         i,q,r+r_offset,s+s_offset,*tmp);
                  tmp++;
                  }
                }
              }
            }
          lock->unlock();
          }
#endif

        timer->enter("2. q.t.");
        // Begin second quarter transformation;
        // generate (iq|jr) for i active and j active or frozen
        for (i=0; i<ni; i++) {
          for (j=0; j<nocc; j++) {

            bzerofast(iqjs_contrib, nbasis*nfuncmax);
            bzerofast(iqjr_contrib, nbasis*nfuncmax);

            for (bf1=0; bf1<ns; bf1++) {
              s = s_offset + bf1;
              double *c_sj = &scf_vector[s][j];
              double *iqjr_ptr = iqjr_contrib;
              for (bf2=0; bf2<nr; bf2++) {
                r = r_offset + bf2;
                if (r > s) {
                  break; // skip to next bf1 value
                  }
                double c_rj = scf_vector[r][j];
                iqjs_ptr = &iqjs_contrib[bf1*nbasis];
                iqrs_ptr = &integral_iqrs[bf1 + ns*nbasis*(bf2 + nr*i)];
                for (q=0; q<nbasis; q++) {
                  *iqjs_ptr++ += c_rj * *iqrs_ptr;
                  if (r != s) *iqjr_ptr += *c_sj * *iqrs_ptr;
                  iqjr_ptr++;
                  iqrs_ptr += ns;
                  } // exit q loop
                // every so often process outstanding messages
                if ((catchup_ctr++ & catchup_mask) == 0) { lock->lock(); mem->catchup(); lock->unlock(); }
                }   // exit bf2 loop
              }     // exit bf1 loop

            // We now have contributions to iqjs and iqjr for one pair i,j,
            // all q, r in R and s in S; send iqjs and iqjr to the node
            // (ij_proc) which is going to have this ij pair
            int ij_proc =  (i*nocc + j)%nproc;
            int ij_index = (i*nocc + j)/nproc;

            // Sum the iqjs_contrib to the appropriate place
            int ij_offset = nbasis*(s_offset + nbasis*ij_index);
            lock->lock();
            mem->sum_reduction_on_node(iqjs_contrib,
                                       ij_offset, ns*nbasis, ij_proc);

            ij_offset = nbasis*(r_offset + nbasis*ij_index);
            mem->sum_reduction_on_node(iqjr_contrib,
                                       ij_offset, nr*nbasis, ij_proc);
            lock->unlock();

            }     // exit j loop
          }       // exit i loop
        // end of second quarter transformation
        timer->exit("2. q.t.");

        }     // endif
      }       // exit R loop
    }         // exit S loop

  lock->lock();
  delete[] integral_iqrs;
  delete[] iqjs_contrib;
  delete[] iqjr_contrib;
  lock->unlock();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
