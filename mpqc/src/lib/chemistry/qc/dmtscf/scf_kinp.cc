
extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <tmpl.h>

#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>

#include <math/array/math_lib.h>
#include <util/misc/libmisc.h>
}

#include <util/keyval/keyval.h>

#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <scf_dmt.h>


///////////////////////////////////////////////////////////////////////////

static char *
newstring(const char *s)
{
  if (!s) return 0;

  char *ret = new char[strlen(s)+1];
  strcpy(ret,s);
  return ret;
}

///////////////////////////////////////////////////////////////////////////
//
// use this if you're lazy.  given just a keyval, initialize centers,
// sym_info, and scf_info by reading the input
//
// input:
//   keyval = a reference to a keyval (preferably a parsed one)
//   centers = reference to uninitialized centers struct
//   scf_info = reference to uninitialized scf struct
//   sym_info = reference to uninitialized sym struct
//
// on return:
//   centers, scf_info, and sym_info are completely initialized
//
// return 0 on success, -1 on failure
//

int
scf_init_scf(KeyVal& keyval, centers_t& centers, scf_struct_t& scf_info,
             sym_struct_t& sym_info)
{

 // initialize the centers struct and the sym_struct
  if (sym_init_centers(keyval,centers,sym_info) < 0) {
    fprintf(stderr,"scf_init_scf:  trouble in sym_init_centers_kv\n");
    return -1;
  }

 // and then fill in the scf struct
  if (scf_init_scf_struct(keyval,centers,scf_info) < 0) {
    fprintf(stderr,"scf_init_scf:  trouble in scf_init_scf_struct\n");
    return -1;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////
//
// given an initialized centers struct and a keyval, read from the keyval
// to initialize an scf struct
//
// input:
//   keyval = a reference to a keyval (preferably a parsed one)
//   centers = reference to an initialized centers struct
//   scf_info = reference to uninitialized scf struct
//
// on return:
//   scf_info contains all necessary information from the input
//
// return 0 on success, -1 on failure
//

int
scf_init_scf_struct(KeyVal& keyval, centers_t& centers,scf_struct_t& scf_info)
{
  int i;
  double nuclear_charge;

 //
 // first off, initialize the scf struct
 //
  init_scf_struct(&scf_info);

 //
 // now find out if this is open shell or not
 //

  char* open = keyval.pcharvalue("opentype");
  if (keyval.error() != KeyVal::OK) open = newstring("none");

  if (keyval.error() == KeyVal::OK) {
    if (!strcmp(open,"highspin") || /* high-spin open-shell */
        !strcmp(open,"singlet")  || /* open-shell singlet */
        !strcmp(open,"twocon")   || /* TCSCF */
        !strcmp(open,"special")) {  /* alpha and beta coeffs given */

      scf_info.iopen = 1;

     // read in the number of singly occupied MOs.  this is, afterall, an
     // open-shell calculation, so there should be some, right?
      scf_info.nopen = keyval.intvalue("socc");
      if (keyval.error() != KeyVal::OK) {
        fprintf(stderr,"scf_init_scf_struct: "
                       "opentype is %s but there is no socc",open);
        return -1;
      }

      if (!strcmp(open,"highspin")) {
        scf_info.hsos=1;
      } else if (!strcmp(open,"singlet")) {
        scf_info.singlet=1;
      } else if (!strcmp(open,"twocon")) {
        scf_info.twocon=1;
      } else if (!strcmp(open,"special")) {
        scf_info.special=1;
      }
    } else if (strcmp(open,"none")) {
      fprintf(stderr,"scf_init_scf_struct:  unrecognized OPENTYPE: %s\n",open);
      return -1;
    }
  }

  delete[] open;

 //
 // get the directory to be used for the checkpoint files.  use strdup
 // since we don't want to use new for anything contained in the scf struct
 //

  char *ckptdir = keyval.pcharvalue("ckpt_dir");
  if (keyval.error() != KeyVal::OK) ckptdir = newstring("./");
  scf_info.ckptdir = strdup(ckptdir);
  delete[] ckptdir;

 //
 // now get the filename to be used for the checkpoint files
 // this name will be appended to the ckptdir, and will be suffixed
 // .{scfvec,fock,etc...}
 //

  char *fname = keyval.pcharvalue("filename");
  if (keyval.error() != KeyVal::OK) fname = newstring("dmtscf");
  scf_info.fname = strdup(fname);
  delete[] fname;

 //
 // grab the wavefunction type and the derivative level to set the default
 // convergence
 //

  char *wfn=keyval.pcharvalue("wfn");
  if (keyval.error() != KeyVal::OK) wfn = newstring("scf");

  char *dertype=keyval.pcharvalue("dertype");
  if (keyval.error() != KeyVal::OK) dertype = newstring("first");

 // set the convergence.  default to 8 for an energy or a gradient,
 // 10 for non-scf, and 12 for second derivs

  scf_info.convergence = keyval.intvalue("convergence");
  if (keyval.error() != KeyVal::OK) {
    scf_info.convergence = 7;
    if (strcmp(wfn,"scf")) scf_info.convergence = 10;
    if (!strcmp(dertype,"second")) scf_info.convergence = 12;
  }

  delete[] wfn;
  delete[] dertype;

 //
 // if there is an old vector available in a checkpoint file, 
 // use it as an intial guess, if there isn't a converged vector about
 //

  scf_info.restart = keyval.booleanvalue("restart");
  scf_info.warmrestart = keyval.booleanvalue("warmrestart");
  scf_info.proj_vector = keyval.booleanvalue("projected_guess");

 //
 // if the point group is not C1, then perform calculation in the SO basis.
 // this will save time in the diagonalization of the fock matrix, and can
 // lead to greater stability. Currently only subgroups of D2h can be
 // done in the SO basis */
 //
  scf_info.use_symmetry = 0;
#if 0 /* not used currently */
  scf_info.use_symmetry=keyval.booleanvalue("use_symmetry");
#endif

 //
 // more mundane members follow
 //

 // the maximum number of scf iterations
  scf_info.maxiter = keyval.intvalue("maxiter");
  if (keyval.error() != KeyVal::OK) scf_info.maxiter = 40;

 // i don't think these are used anymore
  scf_info.debug = keyval.booleanvalue("debug");
  scf_info.debug_node = keyval.booleanvalue("debug_node");

 // use (nproc-1) nodes for gmat calculation (better load balance)
  scf_info.load_bal = keyval.booleanvalue("load_balance_gmat");

 // add on dft correction?
  scf_info.dft = keyval.booleanvalue("dft");

 // use self-consistent dft? means that we only do J part of Gmat
  scf_info.scdft = keyval.booleanvalue("scdft");

 // should the exchange energy be computed separately?
  scf_info.exchange = keyval.booleanvalue("exchange");

 // eliminate integral batches based on size of pmax?
  scf_info.eliminate = keyval.booleanvalue("eliminate");
  if (keyval.error() != KeyVal::OK) scf_info.eliminate=1;

 // should the checkpoint file be deleted?
  scf_info.ckpt_del = keyval.booleanvalue("ckpt_del");
  if (keyval.error() != KeyVal::OK) scf_info.ckpt_del=1;

 // print flag
  scf_info.print_flg = keyval.intvalue("print_flag");

 // how often to checkpoint
  scf_info.ckpt_freq = keyval.intvalue("ckpt_freq");
  if (keyval.error() != KeyVal::OK) scf_info.ckpt_freq=5;

 // how often to reset the density and fock matrices
  scf_info.p_reset_freq = keyval.intvalue("density_reset_frequency");
  if (keyval.error() != KeyVal::OK) scf_info.p_reset_freq=10;

 // use local density matrices?
  scf_info.local_p = keyval.booleanvalue("local_P");

 //
 // integral elimination a la Alrichs.  I don't yet use minimized density
 // differences, however.
 // if pmax*imax < 10^-int_cut, eliminate that integral batch
 // pmax = MAX( p[ij], p[kl], .25*( p[ik], p[il], p[jk], p[jl])
 // imax = MAX (IJKL), I,J,K,L = shell indices
 //
  scf_info.intcut = keyval.intvalue("threshold");
  if (keyval.error() != KeyVal::OK) scf_info.intcut = 12;

 //
 // this is used by version 2 of the integral library.  set int_store to
 // the number of integrals to be kept in memory.
 //
  scf_info.int_store = keyval.intvalue("integral_storage");

 //
 // the number of error matrices to store in the DIIS procedure.
 // 6 is good for closed-shell, 4 for open-shell
 //
  scf_info.ndiis = keyval.intvalue("ndiis");
  if (keyval.error() != KeyVal::OK)
    scf_info.ndiis = (scf_info.iopen) ? 4 : 6;

 //
 // this is an experimental flag for my own use. it is used to select
 // what to use for the diagonal blocks of the effective fock matrix
 // in an open-shell calculation
 //
  scf_info.fock_typ = keyval.intvalue("fock_type");


 //
 // this is a damping factor for the bmat in the diis procedure.
 // the defaults are pretty darn good
 //
  scf_info.diisdamp = keyval.doublevalue("diisdamp");
  if (keyval.error() != KeyVal::OK) {
    scf_info.diisdamp = (scf_info.iopen) ? 0.02 : 0.0;
    if (scf_info.twocon) scf_info.diisdamp = 0.01;
  }

 // level shifting, of course
  if (scf_info.iopen) {
    scf_info.lvl_shift = keyval.doublevalue("levelshift");
    if (keyval.error() != KeyVal::OK)
      scf_info.lvl_shift = 1.0;
  } else {
    scf_info.lvl_shift = 0.0;
  }

 //
 // find the nuclear charge
 //

  nuclear_charge = 0.0;
  for (i=0; i < centers.n; i++)
    nuclear_charge += centers.center[i].charge;

 // read in the number of mo's to occupy in each symmetry type

  scf_info.nclosed = keyval.intvalue("docc");

  if (keyval.error() != KeyVal::OK && !scf_info.iopen) {
    scf_info.nclosed = (int) ((nuclear_charge+0.5)/2);

    if (fabs(nuclear_charge-2*scf_info.nclosed)> 1.e-5) {
      fprintf(stderr,"Yikes! Error assigning charge--give docc.!\n");
      return -1;
    }
  }

 // use diis extrapolation?
  scf_info.diis_flg = keyval.booleanvalue("diis");
  if (keyval.error() != KeyVal::OK) scf_info.diis_flg = 1;
 
 // what iteration to begin diis extrapolation at
  scf_info.it_diis = keyval.intvalue("diisstart");

 //
 // set up alpha and beta arrays.  These really aren't needed yet.
 // actually, they'll probably never be needed
 //

  if (scf_info.iopen) {
    scf_info.alpha = keyval.doublevalue("alpha");
    if (keyval.error() != KeyVal::OK) scf_info.alpha = 0.0;

    scf_info.beta = keyval.doublevalue("beta");
    if (keyval.error() != KeyVal::OK) scf_info.beta = -1.0;

    scf_info.beta = -scf_info.beta;
  }

  int_initialize_1e(0,0,&centers,&centers);
  int_initialize_offsets1(&centers,&centers);

  int nbasis = centers.nfunc;
  int ntri = nbasis*(nbasis+1)/2;

  scf_info.nbfao = nbasis;
  scf_info.nbatri = ntri;

 // for now, use cartesian functions, so nbfso=nbfao
  scf_info.nbfso = nbasis;
  scf_info.nbstri = ntri;

 //
 // all of this is anachronistic
 // nsomax is the size of the largest symmetry block
 // mxcoef is the sum of the sizes of the symmetry blocks squared
 // mxcoef2 is the ioff version of mxcoef
 //

  scf_info.nsomax = nbasis;
  scf_info.mxcoef = nbasis*nbasis;
  scf_info.mxcoef2 = ntri;

 // deallocate integral stuff
  int_done_offsets1(&centers,&centers);
  int_done_1e();

  return 0;
}
  

/////////////////////////////////////////////////////////////////////////

#define PBOOL(a) ((a) ? ("yes"):("no"))

void
scf_print_options(FILE* outfile, scf_struct_t& scf_info)
{

  fprintf(outfile,"\n\n  math/dmt/libdmtscf options:\n");

  if (!scf_info.iopen) {
    fprintf(outfile,"    opentype         = none\n");
  } else if (scf_info.hsos) {
    fprintf(outfile,"    opentype         = highspin\n");
  } else if (scf_info.singlet) {
    fprintf(outfile,"    opentype         = singlet\n");
  } else if (scf_info.twocon) {
    fprintf(outfile,"    opentype         = twocon\n");
  } else if (scf_info.special) {
    fprintf(outfile,"    opentype         = special\n");
  }

  fprintf(outfile,"    use_symmetry     = %s\n",PBOOL(scf_info.use_symmetry));
  fprintf(outfile,"    restart          = %s\n",PBOOL(scf_info.restart));
  fprintf(outfile,"    warmrestart      = %s\n",PBOOL(scf_info.warmrestart));
  fprintf(outfile,"    projected_guess  = %s\n",PBOOL(scf_info.proj_vector));
  fprintf(outfile,"    diis             = %s\n",PBOOL(scf_info.diis_flg));
  fprintf(outfile,"    ckpt_del         = %s\n",PBOOL(scf_info.ckpt_del));
  fprintf(outfile,"    local_P          = %s\n",PBOOL(scf_info.local_p));
  fprintf(outfile,"    eliminate        = %s\n",PBOOL(scf_info.eliminate));
  fprintf(outfile,"    exchange         = %s\n",PBOOL(scf_info.exchange));
  fprintf(outfile,"    load_balance_gmat= %s\n",PBOOL(scf_info.load_bal));
  fprintf(outfile,"    scdft            = %s\n",PBOOL(scf_info.scdft));
  fprintf(outfile,"    dft              = %s\n",PBOOL(scf_info.dft));
  fprintf(outfile,"    debug            = %s\n",PBOOL(scf_info.debug));
  fprintf(outfile,"    debug_node       = %s\n",PBOOL(scf_info.debug_node));

  fprintf(outfile,"    ckpt_dir         = %s\n",scf_info.ckptdir);
  fprintf(outfile,"    filename         = %s\n",scf_info.fname);

  fprintf(outfile,"    ckpt_freq        = %d\n",scf_info.ckpt_freq);
  fprintf(outfile,"    convergence      = %d\n",scf_info.convergence);
  fprintf(outfile,"    maxiter          = %d\n",scf_info.maxiter);
  fprintf(outfile,"    threshold        = %d\n",scf_info.intcut);
  fprintf(outfile,"    integral_storage = %d\n",scf_info.int_store);
  fprintf(outfile,"    ndiis            = %d\n",scf_info.ndiis);
  fprintf(outfile,"    diisstart        = %d\n",scf_info.it_diis);
  fprintf(outfile,"    print_flag       = %d\n",scf_info.print_flg);
  fprintf(outfile,"    docc             = %d\n",scf_info.nclosed);
  fprintf(outfile,"    socc             = %d\n",scf_info.nopen);

  fprintf(outfile,"    diisdamp         = %f\n",scf_info.diisdamp);
  fprintf(outfile,"    levelshift       = %f\n",scf_info.lvl_shift);

  fprintf(outfile,"    density_reset_frequency = %d\n",scf_info.p_reset_freq);

  fprintf(outfile,"\n  libdmtscf calculated the following\n");
  fprintf(outfile,"    nbasis           = %d\n",scf_info.nbfao);
  fprintf(outfile,"    nbasis (SO)      = %d\n",scf_info.nbfso);
  fprintf(outfile,"    nbatri           = %d\n",scf_info.nbatri);
  fprintf(outfile,"    nbstri           = %d\n",scf_info.nbstri);

  fprintf(outfile,"\n");

  switch (scf_info.fock_typ) {
    case 0:
      break;
    case 1:
      fprintf(outfile,"\n  a fock matrix for high spin will be used\n");
      fprintf(outfile,"  this form may not work well with diis\n");
      break;
    default:
      fprintf(outfile,"\n  an experimental fock matrix will be used\n");
  }

  if (scf_info.iopen) {
    fprintf(outfile,"\n  open-shell energy coeffs\n");
    fprintf(outfile,"  open shell pair    alpha         beta\n");

    fprintf(outfile,"        %d  %d       %f     %f\n",1,1,
                 scf_info.alpha,-scf_info.beta);
  }
}

/////////////////////////////////////////////////////////////////////////

/*****************************************************************************
 *
 * construct a centers struct using the old basis set.  
 *
 * input:
 *   centers = pointer to initialized centers struct
 *   oldcenters = pointer to uninitialized centers struct
 *
 * on return:
 *   oldcenters contains copy of centers but with old basis set
 *
 * return 0 on success, -1 on failure
 */

static void
input_errors(const char* msg)
{
  int errcod = -1;
  bcast0(&errcod,sizeof(int),mtype_get(),0);
  if (msg) fprintf(stderr,"%s\n",msg);
}

int
scf_make_old_centers(KeyVal& topkeyval, centers_t& centers,
                                        centers_t& oldcenters)
{
  int i;
  int errcod;

  if (mynode0() == 0) {
   // create a keyval which will look in the :project section of the input
    PrefixKeyVal pkv(":project",topkeyval); pkv.unmanage();
    AggregateKeyVal keyval(pkv,topkeyval); keyval.unmanage();

   // read the value of oldbasis.
    char *oldbasis = keyval.pcharvalue("oldbasis");
    if (keyval.error() != KeyVal::OK) {
      input_errors("scf_make_old_centers: there is no oldbasis");
      return -1;
    }

   // now allocate memory for oldcenters
    errcod=allocbn_centers(&oldcenters,"n",centers.n);
    if (errcod != 0) {
      input_errors("scf_make_old_centers: could not allocate oldcenters");
      return -1;
    }

  // and then copy much of centers into oldcenters, but hold off on the
  // basis set info
    for (i=0; i < centers.n; i++) {
      center_t *center = &oldcenters.center[i];

      errcod = allocbn_center(center,"atom charge",
                              centers.center[i].atom,
                              centers.center[i].charge);
      if (errcod!=0) {
        fprintf(stderr,"scf_make_old_centers: could not alloc center %d\n",i);
        input_errors(0);
        return -1;
      }

      center->r[0]=centers.center[i].r[0];
      center->r[1]=centers.center[i].r[1];
      center->r[2]=centers.center[i].r[2];

      if (int_read_basis(keyval,sym_to_atom(center->atom),
                                        oldbasis,center->basis) < 0) {
        fprintf(stderr,"scf_make_old_centers: could not read basis %d\n",i);
        input_errors(0);
        return -1;
      }
    }

    int_normalize_centers(&oldcenters);

    errcod=0;
    bcast0(&errcod,sizeof(int),mtype_get(),0);

  /* broadcast old centers struct to nodes */

    bcast0_centers(&oldcenters,0,0);

   /* if we get here all is done */

    return 0;

 /* the other nodes just sit around and wait to see what happened */
  } else {
    bcast0(&errcod,sizeof(int),mtype_get(),0);
  }

  if (errcod != 0) return -1;

 /* get oldcenters from node 0 */
  bcast0_centers(&oldcenters,0,0);

  return 0;
}
