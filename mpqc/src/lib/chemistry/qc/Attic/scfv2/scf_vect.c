#define TIME 0

/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/04/17  15:02:15  seidl
 * add vector projection stuff
 *
 * Revision 1.2  1992/04/06  12:50:08  seidl
 * fix MO message
 *
 * Revision 1.1.1.1  1992/03/17  17:09:22  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:20  seidl
 * Initial revision
 *
 * Revision 1.10  1992/01/16  19:49:40  seidl
 * use new integral routines, call int_done_x when cleaning up
 *
 * Revision 1.9  1992/01/13  19:17:07  seidl
 * add timing statements, replace scf_make_g_d and scf_make_pk with
 * scf_make_gmat
 *
 * Revision 1.8  1992/01/09  14:42:25  seidl
 * print statement causes icc to choke
 *
 * Revision 1.7  1992/01/09  11:53:10  seidl
 * add parallel code
 *
 * Revision 1.6  1992/01/02  18:12:30  seidl
 * have scf_file place total scf energy from the last calculation in
 * the current scf_struct
 *
 * Revision 1.5  1992/01/02  16:23:28  seidl
 * a few cosmetic changes
 *
 * Revision 1.4  1991/12/24  19:34:03  seidl
 * use new zero_ routines, get rid of zerov and zerom
 *
 * Revision 1.3  1991/12/24  11:49:41  seidl
 * add pretty print stuff
 *
 * Revision 1.2  1991/12/20  16:32:19  seidl
 * finish initial revision
 *
 * Revision 1.1  1991/12/17  21:45:25  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/02  19:58:51  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/bio/libbio.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"
#include "scfallc.h"
#include "scfbrd.h"
#include "scfbwr.h"
#include "scffree.h"
#include "scfinit.h"
#include "scfprnt.h"

#include "scf_oeis.gbl"
#include "scf_mkden.gbl"
#include "scf_core.gbl"
#include "scf_iter.gbl"
#include "scf_gmat.gbl"

#include "scf_vect.gbl"
#include "scf_vect.lcl"

GLOBAL_FUNCTION int
scf_vector(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,_cls,_ops,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
scf_mf_ptrs_t *_mfp;
int _mast;
int _cls;
int _ops;
FILE *_outfile;
{
  int i,j,k,l;
  int *mptr=&_mfp->cur_ptr;
  int old_vec;
  int errcod;
  int nbasis=_scf_info->nbfao;
  int ntri=_scf_info->nbatri;

  sym_d_matrix_t scf_vec;
  sym_d_vector_t occ_num,evals;
  sym_d_vector_t pmat,pmato,dpmat,dpmato;
  double_vector_t gmat,gmato;
  scf_pk_buf_c_t clbuf;
  scf_pk_buf_o_t opbuf;

/* calculate one-electron integrals and place them in the master file */
  errcod = scf_oeis(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble in forming one-electron integrals\n");
    return(-1);
    }

/* check to see if there is an old vector, if so read it in, else form an
 * initial guess from the core hamiltonian 
 */

  errcod = alloc_s_d_matrix(&scf_vec,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"could not allocate memory for scf_vec matrix\n");
    return(-1);
    }

  if(_scf_info->restart) {
    scf_struct_t tmp;

/* is there a vector in the house? */
    bio_rewind(_mast);
    errcod = bio_read(_mast,&old_vec,sizeof(int));
    if(old_vec==0) _scf_info->restart=0;
    else {
      *mptr = _mfp->scf_vector;
      free_sym_d_matrix(&scf_vec);
      errcod = bread_sym_d_matrix(_mast,&scf_vec,mptr);
      if(errcod < 0) {
        fprintf(_outfile,"scf_vector:\n");
        fprintf(_outfile,"trouble reading scf vector from master file\n");
        return(-1);
        }
      init_scf_struct(&tmp);
      *mptr = _mfp->scf_struct;
      errcod = bread_scf_struct(_mast,&tmp,mptr);
      if(errcod < 0) {
        fprintf(_outfile,"scf_vector:\n");
        fprintf(_outfile,"trouble reading scf struct from master file\n");
        return(-1);
        }

/* in scf_file the total scf energy from the last iteration is put 
 * into e_elec */
      fprintf(_outfile,"\n  reading old vector from master file\n");
      fprintf(_outfile,"  scf energy from old vector = %20.10f\n\n",
                tmp.e_elec);
      }
    }

  if(!_scf_info->restart) {
    fprintf(_outfile,
      "\n  first run, so defaulting to core-hamiltonian guess\n\n");
    errcod = scf_core_guess(&scf_vec,_irreps,_scf_info,_mfp,_mast,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_vector:\n");
      fprintf(_outfile,"trouble forming guess scf vector\n");
      return(-1);
      }
    if(_scf_info->proj_guess) {
      errcod=project_old_vector(_centers,_scf_info,_sym_info,_irreps,_mfp,_mast,
               &scf_vec,_outfile);
      }
    }

  errcod = alloc_s_d_vector(&pmat,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"could not allocate memory for density matrix\n");
    return(-1);
    }
  errcod = alloc_s_d_vector(&dpmat,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,
      "could not allocate memory for density difference matrix\n");
    return(-1);
    }


  if(_scf_info->iopen) {
    errcod = alloc_s_d_vector(&pmato,_irreps);
    if(errcod != 0) {
      fprintf(_outfile,"scf_vector:\n");
      fprintf(_outfile,
        "could not allocate memory for open-shell density matrix\n");
      return(-1);
      }
    errcod = alloc_s_d_vector(&dpmato,_irreps);
    if(errcod != 0) {
      fprintf(_outfile,"scf_vector:\n");
      fprintf(_outfile,
       "could not allocate memory for open-shell density difference matrix\n");
      return(-1);
      }
    }
  else {
    init_sym_d_vector(&pmato);
    init_sym_d_vector(&dpmato);
    }

/* set up occupation numbers and initialize eigenvalues */

  init_sym_d_vector(&occ_num);
  init_sym_d_vector(&evals);
  occ_num.nirrep=_irreps->nirrep;
  occ_num.ir = 
    (double_vector_t *) malloc(sizeof(double_vector_t)*occ_num.nirrep);
  if(occ_num.ir==NULL) {
    fprintf(_outfile,"could not allocate memory for occ_num vector\n");
    return(-1);
    }

  evals.nirrep=_irreps->nirrep;
  evals.ir = 
    (double_vector_t *) malloc(sizeof(double_vector_t)*evals.nirrep);
  if(evals.ir==NULL) {
    fprintf(_outfile,"could not allocate memory for evals vector\n");
    return(-1);
    }

  for(i=0; i < occ_num.nirrep ; i++) {
    errcod = allocbn_double_vector(&occ_num.ir[i],"n",_irreps->ir[i].num_so);
    if(errcod!=0) {
      fprintf(_outfile,"could not allocate memory for occ_num vector %d\n",i);
      return(-1);
      }
    if(_irreps->ir[i].nclosed+_irreps->ir[i].nopen > _irreps->ir[i].num_so) {
      fprintf(_outfile,"too many electrons in irrep %d(%s)\n",
              i,_irreps->ir[i].irrep_label);
      return(-1);
      }
    for(j=0; j < _irreps->ir[i].nclosed ; j++) occ_num.ir[i].d[j]=2.0;
    for(; j < _irreps->ir[i].nclosed+_irreps->ir[i].nopen ; j++) 
                                               occ_num.ir[i].d[j]=1.0;
    for(; j < _irreps->ir[i].num_so ; j++) occ_num.ir[i].d[j]=0.0;

    errcod = allocbn_double_vector(&evals.ir[i],"n",_irreps->ir[i].num_so);
    if(errcod!=0) {
      fprintf(_outfile,"could not allocate memory for evals vector %d\n",i);
      return(-1);
      }
    }

/* form initial density matrices */

  zero_sym_d_vector(&pmat);
  zero_sym_d_vector(&dpmat);
  if(_scf_info->iopen) {
    zero_sym_d_vector(&pmato);
    zero_sym_d_vector(&dpmato);
    }

  errcod = scf_make_density(_scf_info,_irreps,
      &scf_vec,&pmat,&dpmat,&pmato,&dpmato,&occ_num,_outfile);
  if(errcod!=0) {
    fprintf(_outfile,"trouble forming density matrices\n");
    return(-1);
    }

/* and let's form the pk file */

  errcod = allocbn_double_vector(&gmat,"n",_scf_info->nbatri);
  if(errcod!=0) {
    fprintf(_outfile,"could not allocate memory for gmat %d\n",i);
    return(-1);
    }
  zero_double_vector(&gmat);
  if(_scf_info->iopen) {
    errcod = allocbn_double_vector(&gmato,"n",_scf_info->nbatri);
    if(errcod!=0) {
      fprintf(_outfile,"could not allocate memory for gmat %d\n",i);
      return(-1);
      }
    zero_double_vector(&gmato);
    }

  _scf_info->maxbufc=4096;
  _scf_info->maxbufo=4096;

  errcod = allocbn_scf_pk_buf_c(&clbuf,"maxbufc",_scf_info->maxbufc);
  if(errcod!=0) {
    fprintf(_outfile,"could not allocate memory for clbuf\n");
    return(-1);
    }
  if(_scf_info->iopen) {
    errcod = allocbn_scf_pk_buf_o(&opbuf,"maxbufo",_scf_info->maxbufo);
    if(errcod!=0) {
      fprintf(_outfile,"could not allocate memory for opbuf\n");
      return(-1);
      }
    }

/* form g matrix */
  errcod = scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,
                         _cls,_ops,&gmat,&gmato,&dpmat,&dpmato,&clbuf,&opbuf,
                         1,0,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector: trouble forming gmat\n");
    return(-1);
    }

/* now iterate */

  errcod = scf_iter(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,_cls,_ops,
             &scf_vec,&evals,&occ_num,&gmat,&gmato,&pmat,&dpmat,&pmato,&dpmato,
             &clbuf,&opbuf,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_vector: trouble in scf_iter\n");
    return(-1);
    }

/* pretty print some things to output */

  print_pretty(_outfile,_scf_info,_irreps,_mfp,_mast,
                                           &scf_vec,&evals,&occ_num,&pmat);

/* finally, write scf vector to master file */

  *mptr = _mfp->scf_vector;
  errcod = bwrite_sym_d_matrix(_mast,&scf_vec,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing scf vector to master file\n");
    return(-1);
    }

  _mfp->evals=*mptr;
  errcod = bwrite_sym_d_vector(_mast,&evals,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing evals to master file\n");
    return(-1);
    }

  _mfp->occ_num=*mptr;
  errcod = bwrite_sym_d_vector(_mast,&occ_num,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing occ_num to master file\n");
    return(-1);
    }

  *mptr = _mfp->pmat;
  errcod = bwrite_sym_d_vector(_mast,&pmat,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing scf vector to master file\n");
    return(-1);
    }

  *mptr=_mfp->pmato;
  errcod = bwrite_sym_d_vector(_mast,&pmato,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing scf vector to master file\n");
    return(-1);
    }

  *mptr=_mfp->scf_struct;
  errcod = bwrite_scf_struct(_mast,_scf_info,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing scf struct to master file\n");
    return(-1);
    }

  *mptr=sizeof(int);
  errcod = bwrite_scf_mf_ptrs(_mast,_mfp,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_vector:\n");
    fprintf(_outfile,"trouble writing master file struct to master file\n");
    return(-1);
    }

  bio_rewind(_mast);
  i=1;
  bio_write(_mast,&i,sizeof(int));

  free_double_vector(&gmat);
  free_scf_pk_buf_c(&clbuf);
  free_sym_d_matrix(&scf_vec);
  free_sym_d_vector(&occ_num);
  free_sym_d_vector(&evals);
  free_sym_d_vector(&pmat);
  free_sym_d_vector(&dpmat);
  if(_scf_info->iopen) {
    free_double_vector(&gmato);
    free_scf_pk_buf_o(&opbuf);
    free_sym_d_vector(&pmato);
    free_sym_d_vector(&dpmato);
    }

  int_done_erep();
  int_done_offsets2(_centers,_centers,_centers,_centers);
  int_done_storage();

  return(0);
  }


LOCAL_FUNCTION VOID
print_pretty(_outfile,_scf_info,_irreps,_mfp,_mast,scf_vec,evals,occ_num,pmat)
FILE *_outfile;
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
scf_mf_ptrs_t *_mfp;
int _mast;
sym_d_matrix_t *scf_vec;
sym_d_vector_t *evals;
sym_d_vector_t *occ_num;
sym_d_vector_t *pmat;
{
  int i,j;
  int nn;
  int errcod;
  int stest=1,ttest=1;
  int *mptr=&_mfp->cur_ptr;
  double etot,ekin,ovlp,epot,virial,num_elec;
  sym_d_vector_t smat,tmat;

  init_sym_d_vector(&smat);
  init_sym_d_vector(&tmat);

  *mptr = _mfp->overlap;
  errcod = bread_sym_d_vector(_mast,&smat,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"could not read smat from master file\n");
    stest=0;
    }

  *mptr = _mfp->kinetic;
  errcod = bread_sym_d_vector(_mast,&tmat,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"could not read tmat from master file\n");
    ttest=0;
    }

  fprintf(_outfile,"\n\
    **************************\n\
    *                        *\n\
    *   MOLECULAR ORBITALS   *\n\
    *                        *\n\
    **************************\n");

  ekin=ovlp=num_elec=0.0;
  for(i=0; i < _irreps->nirrep ; i++) {
    if(_irreps->ir[i].num_so) {
      fprintf(_outfile,"\n  irrep %s\n",_irreps->ir[i].irrep_label);
      math_print_dmvv(_outfile,&scf_vec->ir[i],&evals->ir[i],&occ_num->ir[i]);

      if(stest) {
        for(j=0; j < pmat->ir[i].n ; j++)
          ovlp += pmat->ir[i].d[j]*smat.ir[i].d[j];
        }
      if(ttest) {
        for(j=0; j < pmat->ir[i].n ; j++)
          ekin += pmat->ir[i].d[j]*tmat.ir[i].d[j];
        }
      for(j=0; j < occ_num->ir[i].n ; j++) num_elec += occ_num->ir[i].d[j];
      }
    }

  etot = _scf_info->nuc_rep + _scf_info->e_elec;

  if(stest) ovlp /= num_elec;
  if(ttest) {
    epot = etot-ekin;
    virial = epot/etot;
    }
  
  if(!_scf_info->converged)
    fprintf(_outfile,"\n\n%8cWarning!  Calculation has not converged!\n",' ');

  fprintf(_outfile,"\n%8ctotal energy       = %20.12f\n",' ',etot);
  if(ttest) {
    fprintf(_outfile,"%8ckinetic energy     = %20.12f\n",' ',ekin);
    fprintf(_outfile,"%8cpotential energy   = %20.12f\n",' ',epot);
    fprintf(_outfile,"%8cvirial theorem     = %20.12f\n",' ',virial);
    }
  if(stest) fprintf(_outfile,"%8cwavefunction norm  = %20.12f\n",' ',ovlp);
  fflush(_outfile);

  if(smat.ir != NULL) free_sym_d_vector(&smat);
  if(tmat.ir != NULL) free_sym_d_vector(&tmat);
  }

LOCAL_FUNCTION int
project_old_vector(_centers,_scf_info,_sym_info,_irreps,_mfp,_mast,
               _scf_vec,_outfile)
centers_t *_centers;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
scf_mf_ptrs_t *_mfp;
int _mast;
sym_d_matrix_t *_scf_vec;
FILE *_outfile;
{
  int i,j;
  int *mptr=&_mfp->cur_ptr;
  int m,n;
  int old_mast,junk,old_ptr;
  sym_d_matrix_t old_vec;
  double_matrix_t basis_overlap,scr1;
  double_matrix_t ov;
  centers_t old_cen;
  scf_mf_ptrs_t old_ptrs;


  old_mast=bio_init_file(0,"old_master");
  bio_open_file(old_mast,"old_master");

  old_ptr=bio_read(old_mast,&junk,sizeof(int));
  bread_scf_mf_ptrs(old_mast,&old_ptrs,&old_ptr);

  old_ptr=old_ptrs.centers;
  bread_centers(old_mast,&old_cen,&old_ptr);

  old_ptr=old_ptrs.scf_vector;
  bread_sym_d_matrix(old_mast,&old_vec,&old_ptr);

  bio_close_file(old_mast,BIO_KEEP);
  bio_done_file(old_mast);

  int_initialize_1e(0,0,&old_cen,_centers);
  int_initialize_offsets1(&old_cen,_centers);

  m=old_cen.nfunc;
  n=_centers->nfunc;

  allocbn_double_matrix(&basis_overlap,"n1 n2",m,n);
  allocbn_double_matrix(&scr1,"n1 n2",n,m);

  ov=old_vec.ir[0];

  int_overlap(&old_cen,_centers,&basis_overlap);

  math_dmxdm_dm(&basis_overlap,1,&ov,0,&scr1,0,n,m,m,0);

  for(i=0; i < _irreps->ir[0].nclosed ; i++) 
    for(j=0; j < n; j++) 
      _scf_vec->ir[0].d[j][i]=scr1.d[j][i];

  {
    sym_d_vector_t scrsv1;
    double_vector_t scrv1;
    double_matrix_t scr2;

    alloc_s_d_vector(&scrsv1,_irreps);
    allocbn_double_vector(&scrv1,"n",n*(n+1)/2);
    allocbn_double_matrix(&scr2,"n1 n2",n,n);

    scf_schmit(_scf_info,_irreps,_mfp,_mast,_scf_vec,
                        &scrsv1,&scr2,&scrv1,1,_outfile);

    free_double_vector(&scrv1);
    free_sym_d_vector(&scrsv1);
    free_double_matrix(&scr2);
    }

  int_done_offsets1(&old_cen,_centers);
  int_done_1e();

  free_double_matrix(&basis_overlap);
  free_double_matrix(&scr1);
  free_sym_d_matrix(&old_vec);
  free_centers(&old_cen);
  }
