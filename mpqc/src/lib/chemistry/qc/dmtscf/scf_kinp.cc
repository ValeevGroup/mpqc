
/* $Log$
 * Revision 1.1  1994/06/08 01:15:11  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * */

static char rcsid[] = "$Id$";

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
}

#include <scf_dmt.h>
#include <util/keyval/keyval.h>

extern "C" {
#include "scf_init.gbl"
}

#define PBOOL(a) ((a) ? ("yes"):("no"))

static int default_value(KeyVal&keyval,char*keyword);

/* This initializes _scf_info, _sym_info, _irreps, and _centers,
 * given the KeyVal input and the symmetry _unique_centers.
 * There are some memory leaks--the return values of keyval.pcharvalue()
 * are not always deleted.
 */
int
scf_get_keyval_input(KeyVal& keyval,
                     centers_t *_unique_centers,
                     scf_struct_t *_scf_info,
                     sym_struct_t *_sym_info,
                     scf_irreps_t *_irreps,
                     centers_t *_centers,
                     FILE *_outfile)
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

 /* first get the point group symbol of the molecule */

  char *pg = keyval.pcharvalue("symmetry");
  if (keyval.error() != KeyVal::OK) pg = "c1";

 /* find out the opentype of the molecule, 
  * and then initialize the scf_struct */

  char* open = keyval.pcharvalue("opentype");
  if (keyval.error() != KeyVal::OK) open = "none";

  init_scf_struct(_scf_info);

  if(keyval.error()==KeyVal::OK) {
    if(!strcmp(open,"highspin")|| /* high-spin open-shell */
            !strcmp(open,"singlet") || /* open-shell singlet */
            !strcmp(open,"twocon") ||  /* TCSCF */
            !strcmp(open,"special")) { /* alpha and beta coeffs given */

      opensh=1;
      nir=0;

 /* check to see if there is a socc vector.  this is afterall an open-shell
  * calculation, so there should be some, right?
  */
      nir = keyval.count("socc");
      if(keyval.error() != KeyVal::OK) {
        fprintf(_outfile,"\n opentype is %s but there is no ",open);
        fprintf(_outfile,"socc vector\n");
        return(-1);
        }

 /* get the number of irreps with open-shell orbitals */
      for(i=nop=0; i < nir; i++) {
        dumo = keyval.intvalue("socc",i);
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

  _scf_info->ckptdir = keyval.pcharvalue("ckpt_dir");
  if(keyval.error()!=KeyVal::OK) {
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

  _scf_info->fname = keyval.pcharvalue("filename");
  if(keyval.error()!=KeyVal::OK) {
    size=strlen("libscfv3")+1;
    _scf_info->fname = (char *) malloc(size);
    check_alloc(_scf_info->fname,"_scf_info->fname");
    strcpy(_scf_info->fname,"libscfv3");
    }

  _scf_info->debug = keyval.booleanvalue("debug");
  _scf_info->debug_node = keyval.booleanvalue("debug_node");

 /* are the coordinates in angstrom units? */
  angs = keyval.booleanvalue("angstrom");

 /* use (nproc-1) nodes for gmat calculation (better load balance) */
  _scf_info->load_bal=keyval.booleanvalue("load_balance_gmat");

 /* should the exchange energy be computed separately? */
  _scf_info->exchange=keyval.booleanvalue("exchange");

 /* cheat by changing the threshold from iteration to iteration? */
  _scf_info->cheat=keyval.booleanvalue("cheat");

 /* eliminate integral batches based on size of pmax? */
  _scf_info->eliminate = keyval.booleanvalue("eliminate");
  if (keyval.error() != KeyVal::OK) _scf_info->eliminate=1;

 /* should the checkpoint file be deleted? */
  _scf_info->ckpt_del=keyval.booleanvalue("ckpt_del");
  if (keyval.error() != KeyVal::OK) _scf_info->ckpt_del=1;

 /* print flag */
  _scf_info->print_flg=keyval.intvalue("print_flag");

 /* how often to checkpoint */
  _scf_info->ckpt_freq=keyval.intvalue("ckpt_freq");
  if (keyval.error() != KeyVal::OK) _scf_info->ckpt_freq=5;

 /* how often to reset the density and fock matrices */
  _scf_info->p_reset_freq=keyval.intvalue("density_reset_frequency");
  if (keyval.error() != KeyVal::OK) _scf_info->p_reset_freq=10;

 /* if there is an old vector available in a checkpoint file, 
  * use it as an intial guess, if there isn't a converged vector about
  */
  _scf_info->restart=0;
  _scf_info->warmrestart=keyval.booleanvalue("warmrestart");

  _scf_info->proj_vector=keyval.booleanvalue("projected_guess");

 /* use local density matrices? */
  _scf_info->local_p=keyval.booleanvalue("local_P");

 /* if the point group is not C1, then perform calculation in the SO basis.
  * this will save time in the diagonalization of the fock matrix, and can
  * lead to greater stability. Currently only subgroups of D2h can be
  * done in the SO basis */
  _scf_info->use_symmetry=0;
#if 0 /* not used currently */
  _scf_info->use_symmetry=keyval.booleanvalue("use_symmetry");
#endif

 /* integral elimination a la Alrichs.  I don't yet use minimized density
  * differences, however.
  * if pmax*imax < 10^-int_cut, eliminate that integral batch
  * pmax = MAX( p[ij], p[kl], .25*( p[ik], p[il], p[jk], p[jl])
  * imax = MAX (IJKL), I,J,K,L = shell indices
  */
  _scf_info->intcut=keyval.intvalue("threshold");
  if (keyval.error() != KeyVal::OK) _scf_info->intcut=12;

 /* this is used by version 2 of the integral library.  set int_store to
  * the number of integrals to be kept in memory.
  */
  _scf_info->int_store=keyval.intvalue("integral_storage");


 /* the number of error matrices to store in the DIIS procedure.
  * 6 is good for closed-shell, 4 for open-shell, 3 for gvb or TCSCF
  */
  _scf_info->ndiis = default_value(keyval,"ndiis");


 /* this is an experimental flag for my own use. it is used to select
  * what to use for the diagonal blocks of the effective fock matrix
  * in an open-shell calculation
  */
  _scf_info->fock_typ = keyval.intvalue("fock_type");


 /* this is a damping factor for the bmat in the diis procedure.
  * the defaults are pretty darn good
  */
  _scf_info->diisdamp = keyval.doublevalue("diisdamp");
  if (keyval.error() != KeyVal::OK) {
      _scf_info->diisdamp = (opensh) ? 0.02 : 0.0;
      if(_scf_info->twocon) _scf_info->diisdamp = 0.01;
    }

 /* level shifting, of course */
  if(opensh) {
      _scf_info->lvl_shift=keyval.doublevalue("levelshift");
      if (keyval.error() != KeyVal::OK) _scf_info->lvl_shift=(opensh)?1.0:0.0;
    }

  errcod = 
    scf_initialize_given_centers(_outfile,_unique_centers,
                                 _centers,_scf_info,_sym_info,_irreps,pg,angs);
  if(errcod != 0) {
    fprintf(_outfile,"trouble in scf_initialize\n");
    return(-1);
    }

/* find the nuclear charge */

  nuclear_charge = 0.0;
  for (i=0; i<_centers->n; i++) {
    nuclear_charge += _centers->center[i].charge;
    }

/* read in the number of mo's to occupy in each symmetry type */

  nir = _irreps->nirrep;
  count = keyval.count("docc");
  if(_scf_info->iopen) {
    count = keyval.count("socc");
    if(keyval.error() != KeyVal::OK) {
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
  if(keyval.error() != KeyVal::OK && !_scf_info->iopen) {
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
      _irreps->ir[i].nclosed = keyval.intvalue("docc",i);
      n_electron += 2*_irreps->ir[i].nclosed;
      if(_scf_info->iopen) {
        _irreps->ir[i].nopen = keyval.intvalue("socc",i);
        n_electron += _irreps->ir[i].nopen;
        }
      }
    }
  
/* pretty-print some info about the calculation */

  char *wfn=keyval.pcharvalue("wfn");
  if (keyval.error() != KeyVal::OK) wfn="scf";

  char *dertype=keyval.pcharvalue("dertype");
  if (keyval.error() != KeyVal::OK) dertype="first";

  fprintf(_outfile,"\n\n  math/dmt/libdmtscf options:\n");
  fprintf(_outfile,"    wfn              = %s\n",wfn);
  fprintf(_outfile,"    dertype          = %s\n",dertype);
  fprintf(_outfile,"    opentype         = %s\n",open);
  fprintf(_outfile,"    warmrestart      = %s\n",PBOOL(_scf_info->warmrestart));
  fprintf(_outfile,"    integral_storage = %d\n",_scf_info->int_store);
  fprintf(_outfile,"    threshold        = %d\n",_scf_info->intcut);
  fprintf(_outfile,"    eliminate        = %s\n",PBOOL(_scf_info->eliminate));

  _scf_info->convergence = default_value(keyval,"convergence");
  fprintf(_outfile,"    convergence      = %d\n",_scf_info->convergence);

  _scf_info->maxiter = default_value(keyval,"maxiter");
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


  _scf_info->diis_flg=keyval.booleanvalue("diis");
  if (keyval.error() != KeyVal::OK) _scf_info->diis_flg=1;
  if(opensh)
    fprintf(_outfile,"  level shift                      = %f\n",
                                                  _scf_info->lvl_shift);
  if(_scf_info->diis_flg) {
    fprintf(_outfile,"  diis scale factor                = %f\n",
                                                  _scf_info->diisdamp+1.0);
    _scf_info->it_diis = default_value(keyval,"diisstart");
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
    _scf_info->alpha=keyval.doublevalue("alpha");
    if (keyval.error() != KeyVal::OK) _scf_info->alpha=0.0;

    _scf_info->beta=keyval.doublevalue("beta");
    if (keyval.error() != KeyVal::OK) _scf_info->beta=-1.0;

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
  if(keyval.booleanvalue("scf_info")) {
    fprintf(_outfile,"\nscf_struct\n");
    print_scf_struct(_outfile,_scf_info);
    fflush(_outfile);
    }
  
/* now print out sym_struct if you wish*/
  if(keyval.booleanvalue("sym_info")) {
    fprintf(_outfile,"\nsym_struct\n");
    print_sym_struct(_outfile,_sym_info);
    fflush(_outfile);
    }
  
/* now print out irreps struct if you wish*/
  if(keyval.booleanvalue("irrep_info")) {
    fprintf(_outfile,"\nscf_irreps\n");
    print_scf_irreps(_outfile,_irreps);
    fflush(_outfile);
    }
  
/* now print out centers struct if you wish*/
  if(keyval.booleanvalue("centers_info")) {
    fprintf(_outfile,"\ncenters\n");
    print_centers(_outfile,_centers);
    fflush(_outfile);
    }
  fflush(_outfile);

  return 0;
  }

static int
default_value(KeyVal&keyval,char*keyword)
{
  int val;

  if(!strcmp(keyword,"convergence")) {
    char *wfn;
    char *der;

    wfn = keyval.pcharvalue("wfn");
    if (keyval.error() != KeyVal::OK) wfn = "scf";
    der = keyval.pcharvalue("dertype");
    if (keyval.error() != KeyVal::OK) der = "first";

    val = keyval.intvalue(keyword);
    if (keyval.error() != KeyVal::OK) {
        val=7;
        if(strcmp(wfn,"scf")) val = 10;
        if(!strcmp(der,"second")) val = 12;
      }
    }
  else if(!strcmp(keyword,"maxiter")) {
    val = keyval.intvalue(keyword);
    if (keyval.error() != KeyVal::OK) val=40;
    }
  else if(!strcmp(keyword,"diisstart")) {
    val = keyval.intvalue(keyword);
    }
  else if(!strcmp(keyword,"ndiis")) {
    char *ot;

    ot = keyval.pcharvalue("opentype");

    val = keyval.intvalue(keyword);
    if (keyval.error() != KeyVal::OK) {
        val = 6;
        if(!ot && strcmp(ot,"none")) val=4;
        if(!ot && !strcmp(ot,"twocon")) val = 3;
      }
    if (ot) delete[] ot;
    }

  return val;
  }
