
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/03/21  00:38:22  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:08  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:06  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.7  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.6  1992/01/13  19:12:22  seidl
 * remove ip stuff for i860 and ncube
 * get rid of node_master
 *
 * Revision 1.5  1992/01/09  11:44:20  seidl
 * add parallel code
 *
 * Revision 1.4  1992/01/02  18:12:30  seidl
 * have scf_file place total scf energy from the last calculation in
 * the current scf_struct
 *
 * Revision 1.3  1991/12/24  19:29:20  seidl
 * do things a little more intelligently
 *
 * Revision 1.2  1991/12/20  16:28:06  seidl
 * change from void function to int
 *
 * Revision 1.1  1991/12/17  21:42:49  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <tmpl.h>
#include <util/ipv2/ip_libv2.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <libbio.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include <sys/stat.h>

#include "scf.h"
#include "scf_file.gbl"
#include "scf_file.lcl"

GLOBAL_FUNCTION int
scf_open_files(si,sym_i,irr,cen,mfp,mf,cl_scr,op_scr,outfile)
scf_struct_t *si;
sym_struct_t *sym_i;
scf_irreps_t *irr;
centers_t *cen;
scf_mf_ptrs_t *mfp;
int *mf;
int *cl_scr;
int *op_scr;
FILE *outfile;
{
  int i;
  int errcod;
  int old_vector;
  int last_rec;
  int *mfptr=&mfp->cur_ptr;

  char master_path[512];
  char master_name[512];
  struct stat master_file;
  scf_irreps_t irrt;
  sym_struct_t smt;
  scf_struct_t scft;
  double_matrix_t scrm;
  sym_d_vector_t scrsv;
  sym_d_matrix_t scrsm;

  init_scf_mf_ptrs(mfp);

/* push the current working keyword list, and add :default and :scf to cwk */

  ip_cwk_push();
  ip_cwk_add(":default");
  ip_cwk_add(":scf");

/* open the master file:
 *    first byte is a flag indicating whether scf_vector is there
 *    then the scf_info struct is written
 *    then scf_vector */

#if 0 /* bag the master file for now */
  get_file_info("master","volume1","%s",master_path);
  get_file_info("master","name","%s",master_name);
  strcat(master_path,master_name);
  strcat(master_path,".master");

  *mf = bio_open_file(0,"master");
  if(*mf < 0) {
    fprintf(outfile,"\n  trouble opening master file\n");
    return(-1);
    }

/* check to see if there is anything in the master file;
 * if empty, set first word to zero and set restart to false.
 * of course, if this is called from a node program, then all
 * of this is moot, so skip it */

  stat(master_path,&master_file);

  if(master_file.st_size==0) {
    int zero=0;
    *mfptr=bio_write(*mf,&zero,sizeof(int));
    bio_rewind(*mf);
    si->restart=0;
    }

  *mfptr = bio_read(*mf,&old_vector,sizeof(int));
  if(old_vector && si->restart) {
    errcod = bread_scf_mf_ptrs(*mf,mfp,mfptr);
    if(errcod < 0) {
      fprintf(stderr,"trouble reading in master file pointers\n");
      return(-1);
      }

  /* see if the same size basis is being used 
     and that the same symmetry is used */

    *mfptr = mfp->sym_struct;
    init_sym_struct(&smt);
    errcod = bread_sym_struct(*mf,&smt,mfptr);
    if(errcod < 0) {
      fprintf(outfile,"trouble reading in old sym_struct\n");
      si->restart=0;
      }

    if(strcmp(smt.point_group,sym_i->point_group)) si->restart=0;
    free_sym_struct(&smt);

    *mfptr = mfp->irreps;
    init_scf_irreps(&irrt);
    errcod = bread_scf_irreps(*mf,&irrt,mfptr);
    if(errcod < 0) {
      fprintf(outfile,"trouble reading in old irreps\n");
      si->restart=0;
      }

    if(!iseq_scf_irreps(&irrt,irr)) si->restart=0;

/* if restart==1 by this point, we probably want to use the old vector,
 *  so let's grab the scf_struct from last iteration, and copy the
 *  total energy into the current scf_struct */
    if(si->restart) {
      init_scf_struct(&scft);
      *mfptr = mfp->scf_struct;
      errcod = bread_scf_struct(*mf,&scft,mfptr);
      if(errcod < 0) {
        fprintf(outfile,"trouble reading in old scf struct\n");
        si->restart=0;
        }
      else si->e_elec=scft.e_elec+scft.nuc_rep;
      }
    }
  else 
    si->restart=0;

  if(!si->restart) {
    int zero=0;
    bio_rewind(*mf);
    *mfptr=bio_write(*mf,&zero,sizeof(int));
    }
    
  *mfptr = sizeof(int);
  errcod = bwrite_scf_mf_ptrs(*mf,mfp,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing scf_mf_ptrs\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->centers = *mfptr;
  errcod = bwrite_centers(*mf,cen,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing centers\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->scf_struct = *mfptr;
  errcod = bwrite_scf_struct(*mf,si,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing scf_struct\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->sym_struct = *mfptr;
  errcod = bwrite_sym_struct(*mf,sym_i,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing sym_struct\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->irreps = *mfptr;
  errcod = bwrite_scf_irreps(*mf,irr,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing irreps\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

/* let's fill up the rest of the master file with junk */
 
  errcod = allocbn_double_matrix(&scrm,"n1 n2",si->nbfao,si->nbfao);
  if(errcod != 0) {
    fprintf(stderr,"scf_file: trouble allocing scrm\n");
    return(-1);
    }
  zero_double_matrix(&scrm);

  mfp->uaoso = *mfptr;
  errcod = bwrite_double_matrix(*mf,&scrm,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing uaoso\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->usoao = *mfptr;
  errcod = bwrite_double_matrix(*mf,&scrm,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing usoao\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }
  free_double_matrix(&scrm);

  errcod = alloc_s_d_vector(&scrsv,irr);
  if(errcod != 0) {
    fprintf(stderr,"scf_file: trouble allocing scrsv\n");
    return(-1);
    }
  zero_sym_d_vector(&scrsv);

  mfp->overlap = *mfptr;
  errcod = bwrite_sym_d_vector(*mf,&scrsv,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing overlap\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->kinetic = *mfptr;
  errcod = bwrite_sym_d_vector(*mf,&scrsv,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing kinetic\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->nuclear = *mfptr;
  errcod = bwrite_sym_d_vector(*mf,&scrsv,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing nuclear\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->hcore = *mfptr;
  errcod = bwrite_sym_d_vector(*mf,&scrsv,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing hcore\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->pmat = *mfptr;
  errcod = bwrite_sym_d_vector(*mf,&scrsv,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing pmat\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  mfp->pmato = *mfptr;
  errcod = bwrite_sym_d_vector(*mf,&scrsv,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing pmato\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  free_sym_d_vector(&scrsv);

  if(!si->restart) {
    errcod = alloc_s_d_matrix(&scrsm,irr);
    if(errcod != 0) {
      fprintf(stderr,"scf_file: trouble allocing scrsv\n");
      return(-1);
      }
    zero_sym_d_matrix(&scrsm);

    mfp->scf_vector = *mfptr;
    errcod = bwrite_sym_d_matrix(*mf,&scrsm,mfptr);
    if(errcod < 0) {
      fprintf(stderr,"trouble writing scf_vector\n");
      fprintf(stderr,"wrote %d bytes\n",-(errcod));
      return(-1);
      }
    free_sym_d_matrix(&scrsm);
    }

 /* write out current file pointers */
  last_rec = *mfptr;
  *mfptr = sizeof(int);
  errcod = bwrite_scf_mf_ptrs(*mf,mfp,mfptr);
  if(errcod < 0) {
    fprintf(stderr,"trouble writing scf_mf_ptrs\n");
    fprintf(stderr,"wrote %d bytes\n",-(errcod));
    return(-1);
    }

  *mfptr = last_rec;

#endif

 /* open scratch files for the supermatrix */

  if (si->write_super) {
    *cl_scr = bio_open_file(0,"cl_scr");
    if(*cl_scr < 0) {
      fprintf(outfile,"\n  trouble opening closed-shell pk-file\n");
      return(-1);
      }
    if(si->iopen) {
      *op_scr = bio_open_file(0,"op_scr");
      if(*op_scr < 0) {
        fprintf(outfile,"\n  trouble opening open-shell pk-file\n");
        return(-1);
        }
      }
    }
    
  ip_cwk_pop();
  return 0;
  }

GLOBAL_FUNCTION VOID
scf_done(si,mf,cl,op)
scf_struct_t *si;
int mf;
int cl;
int op;
{
#if 0
  bio_close_file(mf,BIO_KEEP);
#endif
  if (si->write_super) {
    bio_close_file(cl,BIO_UNLINK);
    if(si->iopen) bio_close_file(op,BIO_UNLINK);
    }
  }
