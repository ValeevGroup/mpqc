
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.3  1992/04/17  15:01:44  seidl
 * add vector projection stuff
 *
 * Revision 1.2  1992/04/06  12:48:56  seidl
 * include stdio.h
 *
 * Revision 1.1.1.1  1992/03/17  17:08:49  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:48  seidl
 * Initial revision
 *
 * Revision 1.5  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.4  1992/01/09  11:44:54  seidl
 * change format of geometry printout
 *
 * Revision 1.3  1991/12/20  18:16:01  seidl
 * the use_symmetry flag now means to use so's, so set it to 0 by default
 *
 * Revision 1.2  1991/12/20  16:30:25  seidl
 * move things around some, change from void to int function
 *
 * Revision 1.1  1991/12/17  21:43:23  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id";

#include <stdio.h>
#include <tmpl.h>
#include <util/ipv2/ip_libv2.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>

#include "scf.h"
#include "scf_init.gbl"

#include "scf_inp.gbl"
#include "scf_inp.lcl"

#define PBOOL(a) ((a) ? ("yes"):("no"))

GLOBAL_FUNCTION int
scf_get_input(si, sym_i, irr, cen, outfile)
scf_struct_t *si;
sym_struct_t *sym_i;
scf_irreps_t *irr;
centers_t *cen;
FILE *outfile;
{
  int i,j,ij;
  int errcod;
  int count;
  int opensh=0;
  int bool;
  int nir,nop,dumo;
  int noptr;

  char *open="none    ";
  char *wfn="scf  ";
  char *dertype="first ";
  char *pg="c1 ";
  char *bdir,*bfile;

  FILE *input;

  sym_d_matrix_t hcore;

 /* push the current working keyword list, then add :default and :scf
  * to cwk list */

  ip_cwk_push();
  ip_cwk_add(":default");
  ip_cwk_add(":scf");

 /* first get the point group symbol of the molecule */

  errcod = ip_string("symmetry",&pg,0);

 /* find out the opentype of the molecule, 
  * and then initialize the scf_struct */

  errcod = ip_string("opentype",&open,0);

  if(errcod==IPE_OK) {
    if(!strcmp(open,"none")) {
      allocbn_scf_struct(si,"optri",0);
      }
    else if(!strcmp(open,"highspin")|| !strcmp(open,"singlet") ||
            !strcmp(open,"twocon") || !strcmp(open,"special")) {

      opensh=1;
      nir=0;
      errcod=ip_count("socc",&nir,0);
      if(errcod!=IPE_OK) {
        fprintf(outfile,"\n opentype is %s but there is no ",open);
        fprintf(outfile,"socc vector\n");
        return(-1);
        }
      for(i=nop=0; i < nir; i++) {
        dumo=0;
        errcod=ip_data("socc","%d",&dumo,1,i);
        if(dumo) nop++;
        }
      noptr=nop*(nop+1)/2;
      if(!strcmp(open,"highspin")) {
        allocbn_scf_struct(si,"hsos n_open optri",
          1,nop,noptr);
        }
      else if(!strcmp(open,"singlet")) {
        allocbn_scf_struct(si,"singlet n_open optri",
          1,nop,noptr);
        }
      else if(!strcmp(open,"twocon")) {
        allocbn_scf_struct(si,"twocon n_open optri",
          1,nop,noptr);
        }
      else if(!strcmp(open,"special")) {
        allocbn_scf_struct(si,"special n_open optri",
          1,nop,noptr);
        }
      }
    else {
      fprintf(outfile," unrecognized OPENTYPE: %s\n",open);
      return(-1);
      }
    }
  else {
    allocbn_scf_struct(si,"optri",0);
    }

  si->proj_guess=0;
  errcod = ip_boolean("project",&si->proj_guess,0);

  si->restart=1;
  errcod = ip_boolean("restart",&si->restart,0);

  si->use_symmetry=0;
  errcod = ip_boolean("use_symmetry",&si->use_symmetry,0);

  si->expense=0;
  errcod= ip_data("expense","%d",&si->expense,0);

  si->intcut=12;
  errcod= ip_data("threshold","%d",&si->intcut,0);

  si->int_store=0;
  errcod= ip_data("integral_storage","%d",&si->int_store,0);

  si->save_thr=5;
  errcod= ip_data("save_thr","%d",&si->save_thr,0);

  si->ndiis = scf_def_value("ndiis");

  si->fock_typ = 0;
  errcod = ip_data("fock_type","%d",&si->fock_typ,0);

  si->diisdamp = (opensh) ? 0.02 : 0.0;
  if(si->twocon) si->diisdamp = 0.01;
  errcod = ip_data("diisdamp","%lf",&si->diisdamp,0);

  si->lvl_shift= (opensh) ? 1.0 : 0.0;
  if(opensh) errcod = ip_data("levelshift","%lf",&si->lvl_shift,0);

/******************************************************************/
 /* open basis set library */

 /* See if the user would like to specify nondefault basis set data bases. */

  ip_cwk_push();
  ip_cwk_add(":default");
  ip_cwk_add(":scf");

  count=0;
  errcod = ip_count("basisfiles",&count,0);
  if(count) ip_append_from_input("basis",outfile);
  else {
    fprintf(outfile," you must supply \"basisdir\" and \"basisfiles\"\n");
    return(-1);
    }

/* read in the centers information from the input */

  errcod = scf_initialize(outfile,cen,si,sym_i,irr,pg);
  if(errcod != 0) {
    fprintf(outfile,"trouble in scf_initialize\n");
    return(-1);
    }

  ip_cwk_pop();

/* read in the number of mo's to occupy in each symmetry type */

  nir = irr->nirrep;
  errcod = ip_count("docc",&count,0);
  if(errcod != 0 && !si->iopen) {
    fprintf(outfile,"Hey! Add some electrons buddy!\n");
    return(-1);
    }
  else {
    if(count != nir) {
      fprintf(outfile,"DOCC vector is the wrong size\n");
      fprintf(outfile,"is %d, should be %d\n",count,nir);
      return(-1);
      }
    }
  if(si->iopen) {
    errcod = ip_count("socc",&count,0);
    if(errcod != 0) {
      fprintf(outfile,"Hey! Add some unpaired electrons buddy!\n");
      return(-1);
      }
    else {
      if(count != nir) {
        fprintf(outfile,"SOCC vector is the wrong size\n");
        fprintf(outfile,"is %d, should be %d\n",count,nir);
        return(-1);
        }
      }
    }

  for(i=0; i < nir ; i++) {
    errcod = ip_data("docc","%d",&irr->ir[i].nclosed,1,i);
    if(si->iopen) {
      errcod = ip_data("socc","%d",&irr->ir[i].nopen,1,i);
      }
    }
  
/* pretty-print some info about the calculation */

  errcod = ip_string("wfn",&wfn,0);
  errcod = ip_string("dertype",&dertype,0);

  fprintf(outfile,"\n\n  wfn         = %s\n",wfn);
  fprintf(outfile,"  dertype     = %s\n",dertype);
  fprintf(outfile,"  opentype    = %s\n",open);
  fprintf(outfile,"  restart     = %s\n",PBOOL(si->restart));
  fprintf(outfile,"  expense     = %d\n",si->expense);
  fprintf(outfile,"  threshold   = %d\n",si->intcut);
  fprintf(outfile,"  save_thr    = %d\n",si->save_thr);

  si->convergence = scf_def_value("convergence");
  fprintf(outfile,"  convergence = %d\n",si->convergence);

  si->maxiter = scf_def_value("maxiter");
  fprintf(outfile,"  maxiter     = %d\n",si->maxiter);
  fprintf(outfile,"  symmetry    = %s\n",pg);
  fprintf(outfile,"  nbasis      = %d\n\n",si->nbfao);

  si->diis_flg=1;
  errcod = ip_boolean("diis",&si->diis_flg,0);
  if(opensh)
    fprintf(outfile,"  level shift                      = %f\n",
                                                  si->lvl_shift);
  if(si->diis_flg) {
    fprintf(outfile,"  diis scale factor                = %f\n",
                                                  si->diisdamp+1.0);
    si->it_diis = scf_def_value("diisstart");
    fprintf(outfile,"  iterations before extrapolation  = %d\n",si->it_diis);
    fprintf(outfile,"  %d error matrices will be kept\n",si->ndiis);
    }
  else fprintf(outfile,"\n  diis turned off\n");

  switch (si->fock_typ) {
  case 0:
    break;
  case 1:
    fprintf(outfile,"\n  a fock matrix for high spin will be used\n");
    fprintf(outfile,"  this form may not work well with diis\n");
    break;
  default:
    fprintf(outfile,"\n  an experimental fock matrix will be used\n");
    }

/* set up alpha and beta arrays.  These really aren't needed yet */

  if(si->iopen) {
    int mm1=1,mm2=1;
    fprintf(outfile,"\n  open-shell energy coeffs\n");
    fprintf(outfile,"  open shell pair    alpha         beta\n");

    if (si->twocon) {
      if(si->n_open == 2) {
        si->alpha[0] = 0.0;
        si->alpha[1] = 0.0;
        si->alpha[2] = 0.0;
        si->beta[0] = 0.0;
        si->beta[1] = 1.0;
        si->beta[2] = 0.0;
        }
      else {
        fprintf(outfile,"This program cannot handle same symmetry");
        fprintf(outfile," tcscf. Try SCFX\n");
        return(-1);
        }
      }
    else if(si->singlet) {
      if(si->n_open == 2) {
        si->alpha[0] = 0.0;
        si->alpha[1] = 0.0;
        si->alpha[2] = 0.0;
        si->beta[0] = 1.0;
        si->beta[1] = -3.0;
        si->beta[2] = 1.0;
        }
      else {
        fprintf(outfile,"This program cannot handle same symmetry");
        fprintf(outfile," singlets. Try SCFX\n");
        return(-1);
        }
      }
    else if(si->hsos) {
      for(i=0; i < si->optri ; i++) {
        si->alpha[i]=0.0;
        si->beta[i]=1.0;
        }
      }
    else {
      for (i=0; i < si->optri ; i++) {
         errcod = ip_data("alpha","%lf",&si->alpha[i],1,i);
         errcod = ip_data("beta","%lf",&si->beta[i],1,i);
         si->beta[i] = -si->beta[i];
         }
      }
    for (i=0; i < si->optri; i++) {
      fprintf(outfile,"        %d  %d       %f     %f\n",mm1,mm2,
                 si->alpha[i],-si->beta[i]);
      mm2++;
      if (mm2 > mm1) {
        mm1++;
        mm2 = 1;
        }
      }
    }

  fprintf(outfile,"\n\n         molecular geometry in atomic units");
  fprintf(outfile,"        basis set\n\n");
  for(i=0; i < cen->n ; i++) 
    fprintf(outfile,"%5d %3s %12.7f %12.7f %12.7f      %s\n",i,
              cen->center[i].atom,
              cen->center[i].r[0],cen->center[i].r[1],cen->center[i].r[2],
              cen->center[i].basis.name);
  fprintf(outfile,"\n");


/**************************************************************/
/* now print out scf_struct if you wish*/
  bool=0;
  errcod = ip_boolean("scf_info",&bool,0);
  if(bool) {
    fprintf(outfile,"\nscf_struct\n");
    print_scf_struct(outfile,si);
    fflush(outfile);
    }
  
/* now print out sym_struct if you wish*/
  bool=0;
  errcod = ip_boolean("sym_info",&bool,0);
  if(bool) {
    fprintf(outfile,"\nsym_struct\n");
    print_sym_struct(outfile,sym_i);
    fflush(outfile);
    }
  
/* now print out irreps struct if you wish*/
  bool=0;
  errcod = ip_boolean("irrep_info",&bool,0);
  if(bool) {
    fprintf(outfile,"\nscf_irreps\n");
    print_scf_irreps(outfile,irr);
    fflush(outfile);
    }
  
/* now print out centers struct if you wish*/
  bool=0;
  errcod = ip_boolean("centers_info",&bool,0);
  if(bool) {
    int_initialize_offsets1(cen,cen);
    fprintf(outfile,"\ncenters\n");
    print_centers(outfile,cen);
    fflush(outfile);
    int_done_offsets1(cen,cen);
    }
  fflush(outfile);

  ip_cwk_pop();

  return 0;
  }

GLOBAL_FUNCTION int
scf_def_value(keyword)
char *keyword;
{
  int val=0;
  int errcod;

  ip_cwk_push();
  ip_cwk_add(":default");
  ip_cwk_add(":scf");

  if(!strcmp(keyword,"convergence")) {
    char *wfn="scf";
    char *der="first";

    errcod = ip_string("wfn",&wfn,0);
    errcod = ip_string("dertype",&der,0);

    val=7;
    if(strcmp(wfn,"scf")) val = 10;
    if(!strcmp(der,"second")) val = 12;
    errcod = ip_data(keyword,"%d",&val,0);
    }
  else if(!strcmp(keyword,"maxiter")) {
    val=40;
    errcod = ip_data(keyword,"%d",&val,0);
    }
  else if(!strcmp(keyword,"diisstart")) {
    errcod = ip_data(keyword,"%d",&val,0);
    }
  else if(!strcmp(keyword,"ndiis")) {
    char *ot="none";

    val = 6;
    errcod = ip_string("opentype",&ot,0);
    if(strcmp(ot,"none")) val=4;
    if(!strcmp(ot,"twocon")) val = 3;
    errcod = ip_data(keyword,"%d",&val,0);
    }

  ip_cwk_pop();
  return val;
  }
