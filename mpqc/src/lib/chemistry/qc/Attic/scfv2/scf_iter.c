#define TIME 0

/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/06/17  22:14:12  jannsen
 * clean up for saber-c
 *
 * Revision 1.2  1992/04/06  12:49:14  seidl
 * include math.h
 *
 * Revision 1.1.1.1  1992/03/17  17:08:52  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:51  seidl
 * Initial revision
 *
 * Revision 1.8  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.7  1992/01/13  19:14:26  seidl
 * add some timing statements, call scf_make_gmat instead of
 * scf_make_g_d and scf_make_g
 *
 * Revision 1.6  1992/01/09  11:58:50  seidl
 * first parallel version
 *
 * Revision 1.5  1992/01/06  11:48:04  seidl
 * add call to scf_schmit
 *
 * Revision 1.4  1992/01/02  16:20:13  seidl
 * add open-shell and diis
 *
 * Revision 1.3  1991/12/24  19:30:41  seidl
 * fix bug in writing of open-shell gmatrix to scr array
 *
 * Revision 1.2  1991/12/24  11:47:47  seidl
 * clean up matrices when converged
 *
 * Revision 1.1  1991/12/20  16:23:08  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/bio/libbio.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include <comm/picl/picl.h>
#include "scf.h"
#include "scfallc.h"
#include "scfbrd.h"
#include "scfbwr.h"
#include "scffree.h"
#include "scfinit.h"
#include "scfprnt.h"

#include "scf_mkden.gbl"
#include "scf_gmat.gbl"
#include "scf_en.gbl"
#include "scf_fock.gbl"
#include "scf_diis.gbl"
#include "scf_orth.gbl"

#include "scf_iter.gbl"
#include "scf_iter.lcl"

GLOBAL_FUNCTION int
scf_iter(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,_cls,_ops,
       _scf_vec,_evals,_occ_num,_gmat,_gmato,_pmat,_dpmat,_pmato,_dpmato,
       _clbuf,_opbuf,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
scf_mf_ptrs_t *_mfp;
int _mast;
int _cls;
int _ops;
sym_d_matrix_t *_scf_vec;
sym_d_vector_t *_evals;
sym_d_vector_t *_occ_num;
double_vector_t *_gmat;
double_vector_t *_gmato;
sym_d_vector_t *_pmat;
sym_d_vector_t *_dpmat;
sym_d_vector_t *_pmato;
sym_d_vector_t *_dpmato;
scf_pk_buf_c_t *_clbuf;
scf_pk_buf_o_t *_opbuf;
FILE *_outfile;
{
  int i,j,k,l;
  int m,nn;
  int iter,maxiter;
  int iopen=_scf_info->iopen;
  int *mptr=&_mfp->cur_ptr;
  int old_vec;
  int errcod;
  int nbasis=_scf_info->nbfao;
  int ntri=_scf_info->nbatri;
  double btime,etime,ttime,clock0();

  double_matrix_t scr1,scr2,scr3;
  double_vector_t scrv1,scrv2;

  sym_d_vector_t fock,focko;
  sym_d_vector_t scrsv1,scrsv2;

  maxiter = _scf_info->maxiter;

/* allocate arrays and such */

  errcod = alloc_s_d_vector(&fock,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: ");
    fprintf(_outfile,"could not allocate memory for fock vector\n");
    return(-1);
    }
  if(iopen) {
    errcod = alloc_s_d_vector(&focko,_irreps);
    if(errcod != 0) {
      fprintf(_outfile,"scf_iter: ");
      fprintf(_outfile,"could not allocate memory for focko vector\n");
      return(-1);
      }
    }
  errcod = alloc_s_d_vector(&scrsv1,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: ");
    fprintf(_outfile,"could not allocate memory for scrsv1 vector\n");
    return(-1);
    }
  errcod = alloc_s_d_vector(&scrsv2,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: ");
    fprintf(_outfile,"could not allocate memory for scrsv2 vector\n");
    return(-1);
    }

  errcod = allocbn_double_matrix(&scr1,"n1 n2",nbasis,nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: trouble allocating scr1\n");
    return(errcod);
    }
  errcod = allocbn_double_matrix(&scr2,"n1 n2",nbasis,nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: trouble allocating scr2\n");
    return(errcod);
    }
  errcod = allocbn_double_vector(&scrv1,"n",ntri);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: trouble allocating scrv1\n");
    return(errcod);
    }
  errcod = allocbn_double_vector(&scrv2,"n",ntri);
  if(errcod != 0) {
    fprintf(_outfile,"scf_iter: trouble allocating scrv2\n");
    return(errcod);
    }

/* begin the iteration here */

#if TIME
  btime = ttime = clock0();
#endif
  for(iter=0; iter < maxiter ; iter++) {

/* form full g matrix from the skeleton gmatrix, place the result in scrv1 */

    if(_sym_info->g > 1) {
#if defined(I860) || defined(NCUBE)
      double_vector_t scrg,scrgo;

      init_double_vector(&scrg);
      init_double_vector(&scrgo);

      assign_double_vector(&scrg,_gmat);
      gsum0(scrg.d,scrg.n,5,33,0);
      bcast0_double_vector(&scrg,0,0);
      sym_sym_vector(_centers,_sym_info,&scrg,&scrv1,&scr1,&scr2,_outfile);
#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for sym_vector %lf\n",etime-btime);
      btime=etime;
#endif

      if(iopen) {
        assign_double_vector(&scrgo,_gmato);
        gsum0(scrgo.d,scrg.n,5,34,0);
        bcast0_double_vector(&scrgo,0,0);
        sym_sym_vector(_centers,_sym_info,&scrgo,&scrv2,&scr1,&scr2,_outfile);
        }
#else
      sym_sym_vector(_centers,_sym_info,_gmat,&scrv1,&scr1,&scr2,_outfile);
      if(iopen)
        sym_sym_vector(_centers,_sym_info,_gmato,&scrv2,&scr1,&scr2,_outfile);
#endif

   /* transform gmat from ao to so basis */

      if(_scf_info->use_symmetry) {
        free_double_matrix(&scr1);
        init_double_matrix(&scr1);
        *mptr = _mfp->uaoso;
        errcod = bread_double_matrix(_mast,&scr1,mptr);
        if(errcod < 0) {
          fprintf(_outfile,"scf_vector:\n");
          fprintf(_outfile,"trouble reading uaoso from master file\n");
          return(-1);
          }

        math_dmxdv_dm(&scr1,1,&scrv1,0,&scr2,0,nbasis,nbasis,nbasis,0);
        math_dmxdm_dv(&scr2,0,&scr1,0,&scrv1,0,nbasis,nbasis,nbasis,0);

        if(iopen) {
          math_dmxdv_dm(&scr1,1,&scrv2,0,&scr2,0,nbasis,nbasis,nbasis,0);
          math_dmxdm_dv(&scr2,0,&scr1,0,&scrv2,0,nbasis,nbasis,nbasis,0);
          }
        }
      }
    else {
#if defined(I860) || defined(NCUBE)
      free_double_vector(&scrv1);
      assign_double_vector(&scrv1,_gmat);
      gsum0(scrv1.d,scrv1.n,5,35,0);
      bcast0_double_vector(&scrv1,0,0);
      if(iopen) {
        free_double_vector(&scrv2);
        assign_double_vector(&scrv2,_gmato);
        gsum0(scrv2.d,scrv2.n,5,36,0);
        bcast0_double_vector(&scrv2,0,0);
        }
#else
      free_double_vector(&scrv1);
      assign_double_vector(&scrv1,_gmat);
      if(iopen) {
        free_double_vector(&scrv2);
        assign_double_vector(&scrv2,_gmato);
        }
#endif
      }
#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for ful gmat %lf\n",etime-btime);
      btime=etime;
#endif

    free_sym_d_vector(&scrsv1);
    init_sym_d_vector(&scrsv1);
    *mptr = _mfp->hcore;
    errcod = bread_sym_d_vector(_mast,&scrsv1,mptr);
    if(errcod < 0) {
      fprintf(_outfile,"scf_iter: ");
      fprintf(_outfile,"trouble reading hcore\n");
      return(errcod);
      }

    make_fock(_scf_info,&scrsv1,&fock,&focko,&scrv1,&scrv2);
    if(iopen) sym_dv_t_svec(&scrv2,&scrsv2);
    scf_calculate_energy(_irreps,_scf_info,_mfp,_mast,&fock,_pmat,_dpmat,
      &scrsv2,_pmato,&scrsv1,iter,_outfile);

#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for fock %lf\n",etime-btime);
      btime=etime;
#endif

/* perform diis extrapolation, returns new fock matrix in ao basis */
    if(_scf_info->diis_flg) {
      errcod = scf_diis(_irreps,_scf_info,&fock,&focko,_scf_vec,_occ_num,iter,
              &scr1,&scr2,&scrv1,&scrv2,_outfile);
      if(errcod != 0) {
        fprintf(_outfile,"scf_iter: trouble in scf_diis\n");
        return(-1);
        }
      }

#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for diis %lf\n",etime-btime);
      btime=etime;
#endif
/* if open-shell, form effective fock matrix used to get new vector */
    if(iopen) {
      errcod = scf_make_fock(_irreps,_scf_info,&fock,&focko,_scf_vec,_occ_num,
                    &scrv1,&scrv2,iter,_outfile);
      if(errcod != 0) {
        fprintf(_outfile,"scf_iter: trouble in scf_make_fock\n");
        return(-1);
        }
      }

/* transform fock matrix to mo basis, then diagonalize to obtain new vector */

    for(m=0; m < _irreps->nirrep ; m++) {
      if(nn=_irreps->ir[m].num_so) {
        if(!iopen) {
          math_dmxdv_dm(&_scf_vec->ir[m],1,&fock.ir[m],0,&scr1,0,nn,nn,nn,0);
          math_dmxdm_dv(&scr1,0,&_scf_vec->ir[m],0,&fock.ir[m],0,nn,nn,nn,0);
          }

        math_diag_dv(&fock.ir[m],&_evals->ir[m],&scr2,1,1.0e-15,1);
        math_dmxdm_dm(&_scf_vec->ir[m],0,&scr2,0,&scr1,0,nn,nn,nn,0);

        if(iopen) {
          double occi;
          double occ0=_occ_num->ir[m].d[0];
       
          for(i=0; i < nn; i++) {
            occi=_occ_num->ir[m].d[i];
            if(occi==occ0 && occi) _evals->ir[m].d[i] += _scf_info->lvl_shift;
            else if(occi) _evals->ir[m].d[i] += 0.5*_scf_info->lvl_shift;
            }
          }

        for(i=0; i < nn; i++)
          for(j=0; j < nn; j++)
            _scf_vec->ir[m].d[i][j] = scr1.d[i][j];
        }
      }
#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for diag %lf\n",etime-btime);
      btime=etime;
#endif

/* if converged, return */

    if(_scf_info->converged) {
      free_sym_d_vector(&fock);
      if(iopen) free_sym_d_vector(&focko);
      free_sym_d_vector(&scrsv1);
      free_sym_d_vector(&scrsv2);
      free_double_matrix(&scr1);
      free_double_matrix(&scr2);
      free_double_vector(&scrv1);
      free_double_vector(&scrv2);
      return(0);
      }

/* orthogonalize new vector */

    errcod = scf_schmit(_scf_info,_irreps,_mfp,_mast,_scf_vec,
                        &scrsv1,&scr1,&scrv1,1,_outfile);
    if(errcod!=0) {
      fprintf(_outfile,"trouble orthogonalizing scf vector\n");
      return(-1);
      }
#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for orth %lf\n",etime-btime);
      btime=etime;
#endif

/* make new density matrices */

    errcod = scf_make_density(_scf_info,_irreps,
      _scf_vec,_pmat,_dpmat,_pmato,_dpmato,_occ_num,_outfile);
    if(errcod!=0) {
      fprintf(_outfile,"trouble forming density matrices\n");
      return(-1);
      }

#if TIME
      etime=clock0();
      if(mynode0()==0)
        fprintf(_outfile,"iter: time for dens %lf\n",etime-btime);
      btime=etime;
#endif

/* and armed with the new density matrix, form new gmatrix */

    errcod = scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,
                           _cls,_ops,_gmat,_gmato,_dpmat,_dpmato,_clbuf,_opbuf,
                           0,iter+1,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_iter: trouble forming gmat\n");
      return(-1);
      }

#if TIME
    etime=clock0();
    if(mynode0()==0)
        fprintf(_outfile,"iter: time for new g %lf\n",etime-btime);
    btime=etime;
#endif

    } /* ad nauseum */
#if TIME
    etime=clock0();
    if(mynode0()==0) fprintf(_outfile,"iter: total time = %lf\n",etime-ttime);
#endif

  return 0;
  }

#define IOFF(a) (a)*((a)+1)/2

LOCAL_FUNCTION VOID
make_fock(_scf_info,_hcore,_fock,_focko,_gmat,_gmato)
scf_struct_t *_scf_info;
sym_d_vector_t *_hcore;
sym_d_vector_t *_fock;
sym_d_vector_t *_focko;
double_vector_t *_gmat;
double_vector_t *_gmato;
{
  int i,j,ij;
  int irrep;
  int iopen=_scf_info->iopen;
  int nso,nsot;
  int joff=0;
  int voff;
  double_vector_t *v,*vo,*h;

 /* F = H + G
  * FO = H + G - GO */

  for(irrep=0; irrep < _fock->nirrep ; irrep++) {
    h = &_hcore->ir[irrep];
    v = &_fock->ir[irrep];
    if(iopen) vo = &_focko->ir[irrep];
    if(nsot=v->n) {
      nso=(int)sqrt(1.0+(double)nsot*8)/2;
      if(iopen) {
        for(i=ij=0; i < nso ; i++) {
          voff=IOFF((joff+i))+joff;
          for(j=0; j <= i ; j++,ij++,voff++) {
            v->d[ij]=h->d[ij]+_gmat->d[voff];
            vo->d[ij] = v->d[ij]-_gmato->d[voff];
            }
          }
        joff += nso;
        }
      else {
        for(i=ij=0; i < nso ; i++) {
          voff=IOFF((joff+i))+joff;
          for(j=0; j <= i ; j++,ij++,voff++) {
            v->d[ij]=h->d[ij]+_gmat->d[voff];
            }
          }
        joff += nso;
        }
      }
    }
  }

LOCAL_FUNCTION VOID
printem(out,vec)
FILE *out;
sym_d_vector_t *vec;
{
  int i;

  for(i=0; i < vec->nirrep ; i++)
    math_print_dv(out,&vec->ir[i]);
  }
