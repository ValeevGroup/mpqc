
/* $Log$
 * Revision 1.2  1994/01/19 13:15:00  seidl
 * add option to use a more load balanced gmat routine.
 *
 * Revision 1.1.1.1  1993/12/29  12:53:15  etseidl
 * SC source tree 0.1
 *
 * Revision 1.15  1992/06/17  21:54:15  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.14  1992/06/16  16:25:51  seidl
 * make p_reset_freq 10
 *
 * Revision 1.13  1992/05/26  20:17:52  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.12  1992/05/19  20:58:27  seidl
 * add cheat stuff
 *
 * Revision 1.11  1992/05/12  10:33:36  seidl
 * no longer print nuclear rep. energy here
 *
 * Revision 1.10  1992/05/04  11:05:11  seidl
 * remove pk ints, add some options
 *
 * Revision 1.9  1992/04/22  15:56:46  seidl
 * add proj_vector
 *
 * Revision 1.8  1992/04/13  11:05:59  seidl
 * alpha and beta no longer necessary
 *
 * Revision 1.7  1992/04/09  17:53:52  seidl
 * indent options
 *
 * Revision 1.6  1992/04/07  18:04:13  jannsen
 * integral_storage is printed out
 *
 * Revision 1.5  1992/04/06  12:36:06  seidl
 * merge in sandia changes
 *
 * Revision 1.4  1992/04/01  01:03:09  seidl
 * fix silly comment end
 *
 * Revision 1.3  1992/03/31  22:24:49  seidl
 * add new options
 *
 * Revision 1.2  1992/03/21  00:42:29  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 * read in boolean eliminate
 *
 * Revision 1.1.1.1  1992/03/17  16:26:20  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:19  seidl
 * Initial revision
 *
 * Revision 1.4  1992/02/26  12:54:38  seidl
 * remember space for null at end of string
 *
 * Revision 1.3  1992/02/18  17:51:44  seidl
 * add local_P option
 *
 * Revision 1.2  1992/02/07  12:59:25  seidl
 * add restart information
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
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

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <tmpl.h>
#include <util/ipv2/ip_libv2.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"
#include "scfinit.h"
#include "scfprnt.h"
#include "scf_init.gbl"

#include "scf_inp.gbl"
#include "scf_inp.lcl"

#define PBOOL(a) ((a) ? ("yes"):("no"))

GLOBAL_FUNCTION int
scf_get_input(_scf_info, _sym_info, _irreps, _centers, _outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
FILE *_outfile;
{
  int i;
  int size;
  int errcod;
  int count;
  int opensh=0;
  int bool;
  int nir,nop,dumo;
  int noptr;
  int angs;
  double nuclear_charge;
  int n_electron;

  char *open="none";
  char *wfn="scf";
  char *dertype="first";
  char *pg="c1";

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

  init_scf_struct(_scf_info);

  if(errcod==IPE_OK) {
    if(!strcmp(open,"highspin")|| /* high-spin open-shell */
            !strcmp(open,"singlet") || /* open-shell singlet */
            !strcmp(open,"twocon") ||  /* TCSCF */
            !strcmp(open,"special")) { /* alpha and beta coeffs given */

      opensh=1;
      nir=0;

 /* check to see if there is a socc vector.  this is afterall an open-shell
  * calculation, so there should be some, right?
  */
      errcod=ip_count("socc",&nir,0);
      if(errcod!=IPE_OK) {
        fprintf(_outfile,"\n opentype is %s but there is no ",open);
        fprintf(_outfile,"socc vector\n");
        return(-1);
        }

 /* get the number of irreps with open-shell orbitals */
      for(i=nop=0; i < nir; i++) {
        dumo=0;
        errcod=ip_data("socc","%d",&dumo,1,i);
        if(dumo) nop++;
        }
      noptr=nop*(nop+1)/2;

      _scf_info->n_open=nop;
      _scf_info->optri=noptr;
      if(!strcmp(open,"highspin")) {
        _scf_info->hsos=1;
        }
      else if(!strcmp(open,"singlet")) {
        _scf_info->singlet=1;
        }
      else if(!strcmp(open,"twocon")) {
        _scf_info->twocon=1;
        }
      else if(!strcmp(open,"special")) {
        _scf_info->special=1;
        }
      }
    else if(strcmp(open,"none")) {
      fprintf(_outfile," unrecognized OPENTYPE: %s\n",open);
      return(-1);
      }
    }


 /* now let's get some other parameters from the input */

 /* get the directory to be used for the checkpoint files. on the i860
  * this should be on the cfs
  */

  errcod = ip_string("ckpt_dir",&_scf_info->ckptdir,0);
  if(errcod!=0) {
#if defined(I860)
    size=strlen("/cfs/")+1;
    _scf_info->ckptdir = (char *) malloc(size);
    check_alloc(_scf_info->ckptdir,"_scf_info->ckptdir");
    strcpy(_scf_info->ckptdir,"/cfs/");
#else
    size=strlen("./")+1;
    _scf_info->ckptdir = (char *) malloc(size);
    check_alloc(_scf_info->ckptdir,"_scf_info->ckptdir");
    strcpy(_scf_info->ckptdir,"./");
#endif
    }

 /* now get the filename to be used for the checkpoint files
  * this name will be appended to the ckptdir, and will be suffixed
  * .{scfvec,fock,etc...}
  */

  errcod = ip_string("filename",&_scf_info->fname,0);
  if(errcod!=0) {
    size=strlen("libscfv3")+1;
    _scf_info->fname = (char *) malloc(size);
    check_alloc(_scf_info->fname,"_scf_info->fname");
    strcpy(_scf_info->fname,"libscfv3");
    }

  _scf_info->debug=0;
  errcod = ip_boolean("debug",&_scf_info->debug,0);
  _scf_info->debug_node=0;
  errcod = ip_boolean("debug_node",&_scf_info->debug_node,0);

 /* are the coordinates in angstrom units? */
  angs=0;
  errcod = ip_boolean("angstrom",&angs,0);

 /* use (nproc-1) nodes for gmat calculation (better load balance) */
  _scf_info->load_bal=0;
  errcod = ip_boolean("load_balance_gmat",&_scf_info->load_bal,0);

 /* should the exchange energy be computed separately? */
  _scf_info->exchange=0;
  errcod = ip_boolean("exchange",&_scf_info->exchange,0);

 /* cheat by changing the threshold from iteration to iteration? */
  _scf_info->cheat=0;
  errcod = ip_boolean("cheat",&_scf_info->cheat,0);

 /* eliminate integral batches based on size of pmax? */
  _scf_info->eliminate=1;
  errcod = ip_boolean("eliminate",&_scf_info->eliminate,0);

 /* should the checkpoint file be deleted? */
  _scf_info->ckpt_del=1;
  errcod = ip_boolean("ckpt_del",&_scf_info->ckpt_del,0);

 /* print flag */
  _scf_info->print_flg=0;
  errcod = ip_data("print_flag","%d",&_scf_info->print_flg,0);

 /* how often to checkpoint */
  _scf_info->ckpt_freq=5;
  errcod = ip_data("ckpt_freq","%d",&_scf_info->ckpt_freq,0);

 /* how often to reset the density and fock matrices */
  _scf_info->p_reset_freq=10;
  errcod = ip_data("density_reset_frequency","%d",&_scf_info->p_reset_freq,0);

 /* if there is an old vector available in a checkpoint file, 
  * use it as an intial guess, if there isn't a converged vector about
  */
  _scf_info->restart=0;
  _scf_info->warmrestart=0;
  errcod = ip_boolean("warmrestart",&_scf_info->warmrestart,0);

  _scf_info->proj_vector=0;
  errcod = ip_boolean("projected_guess",&_scf_info->proj_vector,0);

 /* use local density matrices? */
  _scf_info->local_p=0;
  errcod = ip_boolean("local_P",&_scf_info->local_p,0);

 /* if the point group is not C1, then perform calculation in the SO basis.
  * this will save time in the diagonalization of the fock matrix, and can
  * lead to greater stability. Currently only subgroups of D2h can be
  * done in the SO basis */
  _scf_info->use_symmetry=0;
#if 0 /* not used currently */
  errcod = ip_boolean("use_symmetry",&_scf_info->use_symmetry,0);
#endif

 /* integral elimination a la Alrichs.  I don't yet use minimized density
  * differences, however.
  * if pmax*imax < 10^-int_cut, eliminate that integral batch
  * pmax = MAX( p[ij], p[kl], .25*( p[ik], p[il], p[jk], p[jl])
  * imax = MAX (IJKL), I,J,K,L = shell indices
  */
  _scf_info->intcut=12;
  errcod= ip_data("threshold","%d",&_scf_info->intcut,0);

 /* this is used by version 2 of the integral library.  set int_store to
  * the number of integrals to be kept in memory.
  */
  _scf_info->int_store=0;
  errcod= ip_data("integral_storage","%d",&_scf_info->int_store,0);


 /* the number of error matrices to store in the DIIS procedure.
  * 6 is good for closed-shell, 4 for open-shell, 3 for gvb or TCSCF
  */
  _scf_info->ndiis = scf_def_value("ndiis");


 /* this is an experimental flag for my own use. it is used to select
  * what to use for the diagonal blocks of the effective fock matrix
  * in an open-shell calculation
  */
  _scf_info->fock_typ = 0;
  errcod = ip_data("fock_type","%d",&_scf_info->fock_typ,0);


 /* this is a damping factor for the bmat in the diis procedure.
  * the defaults are pretty darn good
  */
  _scf_info->diisdamp = (opensh) ? 0.02 : 0.0;
  if(_scf_info->twocon) _scf_info->diisdamp = 0.01;
  errcod = ip_data("diisdamp","%lf",&_scf_info->diisdamp,0);

 /* level shifting, of course */
  _scf_info->lvl_shift= (opensh) ? 1.0 : 0.0;
  if(opensh) errcod = ip_data("levelshift","%lf",&_scf_info->lvl_shift,0);

/******************************************************************/
 /* open basis set library. user must provide a basisdir, and one or more
  * basis files.
  */

  ip_cwk_push();
  ip_cwk_add(":default");
  ip_cwk_add(":scf");

  count=0;
  errcod = ip_count("basisfiles",&count,0);
  if(count) ip_append_from_input("basis",_outfile);

/* read in the centers information from the input */

  errcod = 
    scf_initialize(_outfile,_centers,_scf_info,_sym_info,_irreps,pg,angs);
  if(errcod != 0) {
    fprintf(_outfile,"trouble in scf_initialize\n");
    return(-1);
    }

  ip_cwk_pop();

/* find the nuclear charge */

  nuclear_charge = 0.0;
  for (i=0; i<_centers->n; i++) {
    nuclear_charge += _centers->center[i].charge;
    }

/* read in the number of mo's to occupy in each symmetry type */

  nir = _irreps->nirrep;
  errcod = ip_count("docc",&count,0);
  if(_scf_info->iopen) {
    errcod = ip_count("socc",&count,0);
    if(errcod != 0) {
      fprintf(_outfile,"Hey! Add some unpaired electrons buddy!\n");
      return(-1);
      }
    else {
      if(count != nir) {
        fprintf(_outfile,"SOCC vector is the wrong size\n");
        fprintf(_outfile,"is %d, should be %d\n",count,nir);
        return(-1);
        }
      }
    }
  if(errcod != 0 && !_scf_info->iopen) {
    if (nir != 1) {
      fprintf(_outfile,"If there is more than one irrep must give docc!\n");
      return(-1);
      }
    _irreps->ir[0].nclosed = (int) ((nuclear_charge+0.5)/2);
    n_electron = 2*_irreps->ir[0].nclosed;
    if (fabs(nuclear_charge-_irreps->ir[0].nclosed*2)> 1.e-5) {
      fprintf(_outfile,"Yikes! Error assigning charge--give docc.!\n");
      return -1;
      }
    }
  else {
    if(count != nir) {
      fprintf(_outfile,"DOCC vector is the wrong size\n");
      fprintf(_outfile,"is %d, should be %d\n",count,nir);
      return(-1);
      }
    n_electron = 0;
    for(i=0; i < nir ; i++) {
      errcod = ip_data("docc","%d",&_irreps->ir[i].nclosed,1,i);
      n_electron += 2*_irreps->ir[i].nclosed;
      if(_scf_info->iopen) {
        errcod = ip_data("socc","%d",&_irreps->ir[i].nopen,1,i);
        n_electron += _irreps->ir[i].nopen;
        }
      }
    }
  
/* pretty-print some info about the calculation */

  errcod = ip_string("wfn",&wfn,0);
  errcod = ip_string("dertype",&dertype,0);

  fprintf(_outfile,"\n\n  math/dmt/libdmtscf options:\n");
  fprintf(_outfile,"    wfn              = %s\n",wfn);
  fprintf(_outfile,"    dertype          = %s\n",dertype);
  fprintf(_outfile,"    opentype         = %s\n",open);
  fprintf(_outfile,"    warmrestart      = %s\n",PBOOL(_scf_info->warmrestart));
  fprintf(_outfile,"    integral_storage = %d\n",_scf_info->int_store);
  fprintf(_outfile,"    threshold        = %d\n",_scf_info->intcut);
  fprintf(_outfile,"    eliminate        = %s\n",PBOOL(_scf_info->eliminate));

  _scf_info->convergence = scf_def_value("convergence");
  fprintf(_outfile,"    convergence      = %d\n",_scf_info->convergence);

  _scf_info->maxiter = scf_def_value("maxiter");
  fprintf(_outfile,"    maxiter          = %d\n",_scf_info->maxiter);
  fprintf(_outfile,"    symmetry         = %s\n",pg);
  fprintf(_outfile,"    local_P          = %s\n",PBOOL(_scf_info->local_p));
  if(_scf_info->print_flg) 
    fprintf(_outfile,"    print_flag       = %d\n",_scf_info->print_flg);
  if(_scf_info->load_bal) 
    fprintf(_outfile,"    load_balance_gmat= %s\n",PBOOL(_scf_info->load_bal));
  if(_scf_info->exchange) 
    fprintf(_outfile,"    exchange         = %s\n",PBOOL(_scf_info->exchange));
  if(_scf_info->cheat) 
    fprintf(_outfile,"    cheat            = %s\n",PBOOL(_scf_info->cheat));
  if(_scf_info->debug) 
    fprintf(_outfile,"    debug            = %s\n",PBOOL(_scf_info->debug));
  if(_scf_info->debug_node) 
    fprintf(_outfile,"    debug_node       = %s\n",
            PBOOL(_scf_info->debug_node));
  fprintf(_outfile,"    nbasis           = %d\n\n",_scf_info->nbfao);


  _scf_info->diis_flg=1;
  errcod = ip_boolean("diis",&_scf_info->diis_flg,0);
  if(opensh)
    fprintf(_outfile,"  level shift                      = %f\n",
                                                  _scf_info->lvl_shift);
  if(_scf_info->diis_flg) {
    fprintf(_outfile,"  diis scale factor                = %f\n",
                                                  _scf_info->diisdamp+1.0);
    _scf_info->it_diis = scf_def_value("diisstart");
    fprintf(_outfile,
       "  iterations before extrapolation  = %d\n",_scf_info->it_diis);
    fprintf(_outfile,"  %d error matrices will be kept\n",_scf_info->ndiis);
    }
  else fprintf(_outfile,"\n  diis turned off\n");


  switch (_scf_info->fock_typ) {
  case 0:
    break;
  case 1:
    fprintf(_outfile,"\n  a fock matrix for high spin will be used\n");
    fprintf(_outfile,"  this form may not work well with diis\n");
    break;
  default:
    fprintf(_outfile,"\n  an experimental fock matrix will be used\n");
    }

/* set up alpha and beta arrays.  These really aren't needed yet.
 * actually, they'll probably never be needed
 */

  if(_scf_info->iopen) {
    _scf_info->alpha=0.0;
    ip_data("alpha","%lf",&_scf_info->alpha,0);

    _scf_info->beta=-1.0;
    ip_data("beta","%lf",&_scf_info->beta,0);

    _scf_info->beta=-_scf_info->beta;

    fprintf(_outfile,"\n  open-shell energy coeffs\n");
    fprintf(_outfile,"  open shell pair    alpha         beta\n");

    fprintf(_outfile,"        %d  %d       %f     %f\n",1,1,
                 _scf_info->alpha,-_scf_info->beta);
    }

#if 0
  fprintf(_outfile,"\n\n         molecular geometry in atomic units");
  fprintf(_outfile,"        basis set\n\n");
  for(i=0; i < _centers->n ; i++) {
    fprintf(_outfile,"%5d %3s %12.7f %12.7f %12.7f      %s\n",i,
              _centers->center[i].atom,
              _centers->center[i].r[0],
              _centers->center[i].r[1],
              _centers->center[i].r[2],
              _centers->center[i].basis.name);
    }
  fprintf(_outfile,"\n");
#endif

  fprintf(_outfile,"\n  number of electrons is %d\n", n_electron);
  fprintf(_outfile,"  overall charge is % 14.6f\n\n",
          nuclear_charge - n_electron);

  fprintf(_outfile,"\n  density matrices will be reset every %d iterations\n",
           _scf_info->p_reset_freq);
  fprintf(_outfile,
    "\n  scf vector and fock matrices will be checkpointed as\n");
  fprintf(_outfile,
      "    %s%s.{scfvec,fock,focko}\n",_scf_info->ckptdir,_scf_info->fname);
  fprintf(_outfile,"  these files will be updated every %d iterations\n",
         _scf_info->ckpt_freq);
  if(_scf_info->ckpt_del)
    fprintf(_outfile,
      "  these files will be deleted if the calculation is successful\n");


/**************************************************************/
/* now print out scf_struct if you wish*/
  bool=0;
  errcod = ip_boolean("scf_info",&bool,0);
  if(bool) {
    fprintf(_outfile,"\nscf_struct\n");
    print_scf_struct(_outfile,_scf_info);
    fflush(_outfile);
    }
  
/* now print out sym_struct if you wish*/
  bool=0;
  errcod = ip_boolean("sym_info",&bool,0);
  if(bool) {
    fprintf(_outfile,"\nsym_struct\n");
    print_sym_struct(_outfile,_sym_info);
    fflush(_outfile);
    }
  
/* now print out irreps struct if you wish*/
  bool=0;
  errcod = ip_boolean("irrep_info",&bool,0);
  if(bool) {
    fprintf(_outfile,"\nscf_irreps\n");
    print_scf_irreps(_outfile,_irreps);
    fflush(_outfile);
    }
  
/* now print out centers struct if you wish*/
  bool=0;
  errcod = ip_boolean("centers_info",&bool,0);
  if(bool) {
    fprintf(_outfile,"\ncenters\n");
    print_centers(_outfile,_centers);
    fflush(_outfile);
    }
  fflush(_outfile);

  ip_cwk_pop();

  return 0;
  }

GLOBAL_FUNCTION int
scf_def_value(keyword)
char *keyword;
{
  int val=0;

  ip_cwk_push();
  ip_cwk_add(":default");
  ip_cwk_add(":scf");

  if(!strcmp(keyword,"convergence")) {
    char *wfn="scf";
    char *der="first";

    ip_string("wfn",&wfn,0);
    ip_string("dertype",&der,0);

    val=7;
    if(strcmp(wfn,"scf")) val = 10;
    if(!strcmp(der,"second")) val = 12;
    ip_data(keyword,"%d",&val,0);
    }
  else if(!strcmp(keyword,"maxiter")) {
    val=40;
    ip_data(keyword,"%d",&val,0);
    }
  else if(!strcmp(keyword,"diisstart")) {
    ip_data(keyword,"%d",&val,0);
    }
  else if(!strcmp(keyword,"ndiis")) {
    char *ot="none";

    val = 6;
    ip_string("opentype",&ot,0);
    if(strcmp(ot,"none")) val=4;
    if(!strcmp(ot,"twocon")) val = 3;
    ip_data(keyword,"%d",&val,0);
    }

  ip_cwk_pop();
  return val;
  }
