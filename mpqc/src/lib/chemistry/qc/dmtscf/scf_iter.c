
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <tmpl.h>
#include <util/group/picl.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_mkden.gbl>
#include <chemistry/qc/dmtscf/scf_gmat.gbl>
#include <chemistry/qc/dmtscf/scf_en.gbl>
#include <chemistry/qc/dmtscf/scf_diis.gbl>
#include <chemistry/qc/dmtscf/scf_orth.gbl>

#include <chemistry/qc/dmtscf/scf_iter.gbl>
#include <chemistry/qc/dmtscf/scf_iter.lcl>

/* these are pointer copies of the core hamiltonian and overlap integrals */
static dmt_matrix Hcore, S;

/* closed shell density and density difference */
static dmt_matrix Pmat, DPmat;

/* open shell density and density difference */
static dmt_matrix PmatO, DPmatO;

/* Gmat and GmatO contain the skeleton G matrices */
static dmt_matrix Gmat, GmatO;

/* scratch matrices...use SScr2 sparingly since it contains the symmetrized
 * open-shell fock matrix after scf_iter() is called
 */
static dmt_matrix Scr1, Scr2, Scr3;
static dmt_matrix SScr1, SScr2;

/**********************************************************************
 *
 * no sense having all these tmp matrices around if you can't use them
 * so here are some functione to grab pointers to them all
 */

GLOBAL_FUNCTION VOID
scf_iter_get_col_tmps(S1, S2, S3)
dmt_matrix *S1;
dmt_matrix *S2;
dmt_matrix *S3;
{
  *S1 = Scr1;
  *S2 = Scr2;
  *S3 = Scr3;
}

GLOBAL_FUNCTION VOID
scf_iter_get_scat_tmps(S1, S2)
dmt_matrix *S1;
dmt_matrix *S2;
{
  *S1 = SScr1;
  *S2 = SScr2;
}

GLOBAL_FUNCTION VOID
scf_iter_get_pmats(S1, S2, S3, S4)
dmt_matrix *S1;
dmt_matrix *S2;
dmt_matrix *S3;
dmt_matrix *S4;
{
  *S1 = Pmat;
  *S2 = DPmat;
  *S3 = PmatO;
  *S4 = DPmatO;
}

GLOBAL_FUNCTION VOID
scf_iter_get_gmats(S1, S2)
dmt_matrix *S1;
dmt_matrix *S2;
{
  *S1 = Gmat;
  *S2 = GmatO;
}


/*************************************************************************
 *
 * given the scf struct this initializes the matrices needed for scf
 * iterations
 *
 * input:
 *   HM = scattered dmt matrix containing core hamiltonian
 *   SM = scattered dmt matrix containing overlap integrals
 */

GLOBAL_FUNCTION int
scf_init_iter(centers,scf_info,HM,SM)
centers_t *centers;
scf_struct_t *scf_info;
dmt_matrix HM;
dmt_matrix SM;
{
  int nbasis = scf_info->nbfao;

  assert(dmt_distribution(HM) == SCATTERED);
  assert(dmt_distribution(SM) == SCATTERED);

  Hcore = HM;
  S = SM;

 /* initialize the diis routines */
  if (scf_init_diis(scf_info) < 0) {
    fprintf(stderr,"scf_init_iter: could not init diis routines\n");
    return -1;
  }

 /* initialize the g matrix stuff */
  if (scf_init_gmat(centers,scf_info) < 0) {
    fprintf(stderr,"scf_init_iter: could not init gmat routines\n");
    return -1;
  }

 /* form initial density matrices */

  Pmat = dmt_create("dmtscf density matrix",nbasis,SCATTERED);
  dmt_fill(Pmat,0.0);

  DPmat = dmt_create("dmtscf density diff matrix",nbasis,SCATTERED);
  dmt_fill(DPmat,0.0);

  if (scf_info->iopen) {
    PmatO = dmt_create("dmtscf open density matrix",nbasis,SCATTERED);
    dmt_fill(PmatO,0.0);

    DPmatO = dmt_create("dmtscf open density diff matrix",nbasis,SCATTERED);
    dmt_fill(DPmatO,0.0);

  } else {
    PmatO = dmt_nil();
    DPmatO = dmt_nil();
  }

 /* allocate memory for the G matrices */

  Gmat = dmt_create("dmtscf g matrix",nbasis,SCATTERED);
  dmt_fill(Gmat,0.0);

  if (scf_info->iopen) {
    GmatO = dmt_create("dmtscf open g matrix",nbasis,SCATTERED);
    dmt_fill(GmatO,0.0);
  } else {
    GmatO = dmt_nil();
  }

 /* allocate scratch arrays now */

  Scr1 = dmt_create("scf_iter: scr1",nbasis,COLUMNS);
  Scr2 = dmt_create("scf_iter: scr2",nbasis,COLUMNS);
  Scr3 = dmt_create("scf_iter: scr3",nbasis,COLUMNS);
  SScr1 = dmt_create("scf_iter: scr4",nbasis,SCATTERED);
  SScr2 = dmt_create("scf_iter: scr5",nbasis,SCATTERED);

  return 0;
}

/************************************************************************
 * 
 * this frees up the matrices used for scf iterations
 *
 * input:
 *   centers  = pointer to initialized centers struct
 *   scf_info = pointer to initialized scf struct
 *   sym_info = pointer to initialized sym struct
 *   Fock     = scattered distributed matrix
 *   FockO    = scattered distributed matrix
 *   Scf_Vec  = column distributed matrix containing the converged wfn
 *
 * on return:
 *   Fock and FockO contain the MO fock matrices
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION VOID
scf_done_iter(centers,scf_info,sym_info,Fock,FockO,Scf_Vec)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
dmt_matrix Fock;
dmt_matrix FockO;
dmt_matrix Scf_Vec;
{

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Fock) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(FockO) == SCATTERED);

 /* clean up the diis functions */
  if (scf_info->diis_flg) scf_done_diis(scf_info);

 /* retrieve the MO Fock matrices */
  make_fock(centers,scf_info,sym_info,Fock,FockO,Gmat,GmatO,SScr1,SScr2);

  dmt_mult(Fock,Scf_Vec,Scr1);
  dmt_mult(Scf_Vec,Scr1,Scr2);
  dmt_copy(Scr2,Fock);

  if (scf_info->iopen) {
    dmt_mult(FockO,Scf_Vec,Scr1);
    dmt_mult(Scf_Vec,Scr1,Scr2);
    dmt_copy(Scr2,FockO);
  }

 /* free up 2ei functions */
  scf_done_gmat(centers,scf_info);

  dmt_free(Pmat);
  dmt_free(DPmat);
  dmt_free(Gmat);
  dmt_free(Scr1);
  dmt_free(Scr2);
  dmt_free(Scr3);
  dmt_free(SScr1);
  dmt_free(SScr2);
  if (scf_info->iopen) {
    dmt_free(PmatO);
    dmt_free(DPmatO);
    dmt_free(GmatO);
  }
}

/*****************************************************************************
 * 
 * actually do the iterations needed for an scf calculation
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   sym_info = pointer to initialized sym struct
 *   centers  = pointer to initialized centers struct
 *   Scf_Vec  = column distributed matrix
 *   HM       = scattered dmt matrix containing core hamiltonian
 *   SM       = scattered dmt matrix containing overlap integrals
 *   Fock     = scattered distributed matrix
 *   FockO    = scattered distributed matrix
 *   evals    = pointer to double array
 *   occ_num  = pointer to double array containing MO occupation numbers
 *   outfile  = FILE pointer to output
 *
 * on return:
 *   Scf_Vec contains the converged SCF wavefunction
 *   Fock and FockO contain the MO fock matrices
 *   evals contains the MO energies
 *
 * return 0 on success, -1 on failure
 */
 
GLOBAL_FUNCTION int
scf_do_iter(scf_info,sym_info,centers,
            Scf_Vec,HM,SM,Fock,FockO,evals,occ_num,outfile)
scf_struct_t *scf_info;
sym_struct_t *sym_info;
centers_t *centers;
dmt_matrix Scf_Vec;
dmt_matrix HM;
dmt_matrix SM;
dmt_matrix Fock;
dmt_matrix FockO;
double *evals;
double *occ_num;
FILE *outfile;
{
  int iter,errcod=0;
  double etot, edif, neelec, delta;
  double plimit = pow(10.0,(double) -(scf_info->convergence));

  char vecfile[512],fockfile[512],fockofile[512];

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(HM) == SCATTERED);
  assert(dmt_distribution(SM) == SCATTERED);
  assert(dmt_distribution(Fock) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(FockO) == SCATTERED);

 /* what are the names of the checkpoint files? */

  sprintf(vecfile,"%s%s.scfvec",scf_info->ckptdir,scf_info->fname);
  sprintf(fockfile,"%s%s.fock",scf_info->ckptdir,scf_info->fname);
  sprintf(fockofile,"%s%s.focko",scf_info->ckptdir,scf_info->fname);

 /* initialize the scf_iter stuff */
  scf_init_iter(centers,scf_info,HM,SM);

 /* print the header for the energies */
  if (mynode0()==0 && outfile) {
    fprintf(outfile,"\n  iter       total energy       "
                    " delta E         delta P          diiser\n");
  }

 /* begin the iteration here */

  scf_info->converged=0;
  scf_info->diis_er=0.0;

  for (iter=0; iter < scf_info->maxiter ; iter++) {

   /* calculate the next scf vector
    * the AO fock matrices are in Fock and FockO now
    */
    errcod = scf_make_ao_fock(scf_info,sym_info,centers,Scf_Vec,Fock,FockO,
                              occ_num,iter,outfile);
    if (errcod < 0) {
      fprintf(stderr,"scf_do_iter:  scf_iter failed on iteration %d\n",iter+1);
      break;
    }

   /* calculate the electronic energy */
    neelec =
      scf_electronic_energy(scf_info,Hcore,Fock,Pmat,FockO,PmatO,SScr1,SScr2);
   
   /* calculate the change in the density */
    delta = scf_rms_delta_density(scf_info,DPmat);

   /* Etot = Eelec + Enucrep */
    etot = scf_info->nuc_rep + neelec;

   /* delta(E) = Eelec(n-1) - Eelec(n) */
    edif = scf_info->e_elec - neelec;
    scf_info->e_elec = neelec;

    if (mynode0()==0 && outfile) {
      fprintf(outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n",
                     iter+1, etot, edif, delta, scf_info->diis_er);
      fflush(outfile);
    }

   /* checkpoint every "ckpt_freq" iterations */
    if (((iter+1) % scf_info->ckpt_freq) == 0) {
      if (mynode0()==0 && outfile) {
        fprintf(outfile,
          "  scf_do_iter: checkpointing vector and fock matrices\n");
      }
      dmt_write(vecfile,Scf_Vec);
      dmt_write(fockfile,Fock);
      if (scf_info->iopen) dmt_write(fockofile,FockO);
    }

   /* if converged, return */
    if (delta < plimit) {
      scf_info->converged = 1;
      if (mynode0()==0 && outfile) {
          fprintf(outfile,"\n  converged scf energy is %20.10f au\n",
                  scf_info->nuc_rep + scf_info->e_elec);
      }
      break;
    }

   /* calculate next vector.  after this, Fock and FockO contain MO matrices */
    errcod = scf_make_new_vec(scf_info,Scf_Vec,Fock,FockO,evals,occ_num,iter);
    if (errcod < 0) {
      fprintf(stderr,"scf_do_iter: could not make new vector %d\n",iter);
      break;
    }
  }

 /* cleanup */

 /* this puts MO fock matrices into Fock and FockO */
  scf_done_iter(centers,scf_info,sym_info,Fock,FockO,Scf_Vec);

  return errcod;
}

/***************************************************************************
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   sym_info = pointer to initialized sym struct
 *   centers  = pointer to initialized centers struct
 *   Scf_Vec  = column distributed matrix containing a trial scf vector
 *   Fock     = scattered distributed matrix
 *   FockO    = scattered distributed matrix
 *   occ_num  = pointer to double array containing MO occupation numbers
 *   iter     = current iteration (first iteration is 0)
 *   outfile  = FILE pointer to output
 *
 * on return:
 *   Fock and FockO contain the AO fock matrices
 *   Fock = Hcore + G
 *   FockO = Hcore + G - GO
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_ao_fock(scf_info,sym_info,centers,Scf_Vec,Fock,FockO,
                                                occ_num,iter,outfile)
scf_struct_t *scf_info;
sym_struct_t *sym_info;
centers_t *centers;
dmt_matrix Scf_Vec;
dmt_matrix Fock;
dmt_matrix FockO;
double *occ_num;
int iter;
FILE *outfile;
{

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Fock) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(FockO) == SCATTERED);

 /* form density matrices from scf vector */
  if (scf_make_density(scf_info,Scf_Vec,Pmat,DPmat,PmatO,DPmatO,occ_num) < 0) {
    fprintf(stderr,"scf_make_ao_fock: trouble forming density matrices\n");
    return -1;
  }

 /* reset density if appropriate */

  if (iter && scf_info->eliminate && ((iter)%scf_info->p_reset_freq == 0)) {
    dmt_copy(Pmat,DPmat);
    dmt_fill(Gmat,0.0);
    if (scf_info->iopen) {
      dmt_copy(PmatO,DPmatO);
      dmt_fill(GmatO,0.0);
    }

    if (mynode0()==0 && outfile)
      fprintf(outfile,"  scf_make_ao_fock: resetting density matrices\n");
  }

 /* and armed with the new density matrix, form new skeleton G and GO */

  if (scf_make_gmat(scf_info,sym_info,centers,
                    Gmat,GmatO,DPmat,DPmatO,SScr1,SScr2,outfile) < 0) {
    fprintf(stderr,"scf_make_ao_fock: trouble forming gmat\n");
    return -1;
  }

  /* form the full G matrices from the skeleton one, and then form 
   * Fock = H + G and FockO = H + G - GO
   */
  make_fock(centers,scf_info,sym_info,Fock,FockO,Gmat,GmatO,SScr1,SScr2);

  return 0;
}

/***************************************************************************
 *
 * input:
 *   scf_info = pointer to initialized scf struct
 *   Scf_Vec  = column distributed matrix
 *   Fock     = scattered distributed matrix
 *   FockO    = scattered distributed matrix
 *   evals    = pointer to double array
 *   occ_num  = pointer to double array containing MO occupation numbers
 *   iter     = current iteration (first iteration is 0)
 *
 * on return:
 *   Scf_Vec contains the current SCF vector
 *   Fock contains the MO fock matrice or effective MO fock matrix (open shell)
 *   FockO contains the open-shell fock matrix transformed to MO basis
 *   evals contains the MO energies
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_make_new_vec(scf_info,Scf_Vec,Fock,FockO,evals,occ_num,iter)
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
dmt_matrix Fock;
dmt_matrix FockO;
double *evals;
double *occ_num;
int iter;
{
  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Fock) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(FockO) == SCATTERED);

 /* perform diis extrapolation, returns new fock matrix in ao basis */
  if (scf_info->diis_flg) {
    if (scf_diis(scf_info,Fock,FockO,Scf_Vec,occ_num,iter,Scr1,Scr2,Scr3) < 0) {
      fprintf(stderr,"scf_make_new_vec: trouble in scf_diis\n");
      return -1;
    }
  }

 /* if open-shell, form effective MO fock matrix used to get new vector
  * and level shift it. otherwise, just transform Fock matrix to MO basis.
  */
  if (scf_info->iopen) {
    make_eff_fock(scf_info,Fock,FockO,Scf_Vec,occ_num,iter);
  } else {
    dmt_mult(Fock,Scf_Vec,Scr1);
    dmt_mult(Scf_Vec,Scr1,Scr2);
    dmt_copy(Scr2,Fock);
  }

 /* now diagonalize the MO Fock matrix */
  dmt_copy(Fock,Scr2); /* dmt_diag needs a columns dist. matrix */
  dmt_diag(Scr2,Scr1,evals);
  dmt_transpose(Scf_Vec);
  dmt_mult(Scf_Vec,Scr1,Scr2);
  dmt_copy(Scr2,Scf_Vec);

 /* un-level shift eigenvalues */

  if (scf_info->iopen) {
    un_level_shift(scf_info,Fock,occ_num,evals);
  }

  if (scf_info->print_flg & 2) {
    sync0();
    tim_print(0);
  }

 /* orthogonalize new vector */
  if (scf_schmidt(scf_info,Scf_Vec,S,1) < 0) {
    fprintf(stderr,"scf_make_new_vec: trouble orthogonalizing scf vector\n");
    return -1;
  }

  return 0;
}

/*************************************************************************
 *
 * this takes the Skeleton fock matrices in Gmat and GmatO, symmetrizes 
 * them, and then forms Fock = Hcore + sym(Gmat), and FockO = sym(GmatO)
 */

LOCAL_FUNCTION VOID
make_fock(centers,scf_info,sym_info,Fock,FockO,Gmat,GmatO,SScr1,SScr2)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
dmt_matrix Fock;
dmt_matrix FockO;
dmt_matrix Gmat;
dmt_matrix GmatO;
dmt_matrix SScr1;
dmt_matrix SScr2;
{
 /* form full g matrix from the skeleton gmatrix, place the result in SScr1 */

  if (sym_info->g > 1) {
    sym_sym_matrix(centers,sym_info,Gmat,SScr1);
    if (scf_info->iopen) sym_sym_matrix(centers,sym_info,GmatO,SScr2);
  } else {
    dmt_copy(Gmat,SScr1);
    if (scf_info->iopen) dmt_copy(GmatO,SScr2);
  }

 /* F = H + G
  * FO = H + G - GO */

  dmt_copy(Hcore,Fock);
  dmt_sum(SScr1,Fock);

  if (scf_info->iopen) {
    dmt_copy(SScr2,FockO);
    dmt_scale(FockO,-1.0);
    dmt_sum(Fock,FockO);
  }
}

LOCAL_FUNCTION VOID
un_level_shift(scf_info,Fock,occ_num,evals)
scf_struct_t *scf_info;
dmt_matrix Fock;
double *occ_num;
double *evals;
{
  int i;
  double *fdiag;
  double occi;
  double occ0 = occ_num[0];

  fdiag = (double *) malloc(sizeof(double)*scf_info->nbfao);
  if (!fdiag) {
    fprintf(stderr,"un_level_shift:  could not malloc fdiag\n");
    exit(1);
  }

  dmt_get_diagonal(Fock,fdiag);
       
  for (i=0; i < scf_info->nbfao; i++) {
    occi = occ_num[i];
    if (occi==occ0 && occi) {
      evals[i] += scf_info->lvl_shift;
      fdiag[i] += scf_info->lvl_shift;
    } else if (occi) {
      evals[i] += 0.5*scf_info->lvl_shift;
      fdiag[i] += 0.5*scf_info->lvl_shift;
    }
  }
  dmt_set_diagonal(Fock,fdiag);
  
  free(fdiag);
}

/*************************************************************************
 * 
 * given the closed and open-shell fock matrices in the AO basis, form
 * the effective open-shell Fock matrix. also level shifts the fock matrix
 *
 * input:
 *   scf_info = pointer to scf struct
 *   Fock     = scattered dmt matrix containing fock matrix (overwritten)
 *   FockO    = scattered dmt matrix containing open-shell fock matrix
 *   Scf_Vec  = column dmt matrix containing scf vector
 *   occ_num  = pointer to double array containing mo occupation numbers
 *   Scr1, Scr2, Scr3 = column dmt scratch matrices
 *   iter     = iteration number
 *
 * on return:
 *   Fock contains effective MO basis fock matrix
 *
 */

LOCAL_FUNCTION VOID
make_eff_fock(scf_info,Fock,FockO,Scf_Vec,occ_num,iter)
scf_struct_t *scf_info;
dmt_matrix Fock;
dmt_matrix FockO;
dmt_matrix Scf_Vec;
double *occ_num;
int iter;
{
  int i,j;
  int ib,isz,ist,jb,jsz,jst;
  int nlocalb=dmt_nlocal(Fock);
  int nn;
  double occi, occj, occ0;
  double *Fblk,*FOblk;

 /* transform fock to mo basis */
  dmt_mult(Fock,Scf_Vec,Scr3);
  dmt_mult(Scf_Vec,Scr3,Scr1);
  dmt_copy(Scr1,Fock);

 /* transform fock_open to mo basis */
  dmt_mult(FockO,Scf_Vec,Scr3);
  dmt_mult(Scf_Vec,Scr3,Scr2);
  dmt_copy(Scr2,FockO);

 /* form effective fock matrix in mo basis */

  occ0 = occ_num[0];
  for (nn=0; nn < nlocalb; nn++ ) {
    dmt_get_block(Fock,nn,&ib,&jb,&Fblk);
    dmt_get_block(FockO,nn,&ib,&jb,&FOblk);
    dmt_describe_block(Fock,ib,&ist,&isz);
    dmt_describe_block(Fock,jb,&jst,&jsz);

    for (i=0; i < isz ; i++) {
      for (j=0; j < jsz ; j++) {
        occi = occ_num[ist+i];
        occj = occ_num[jst+j];

     /* default: Guest & Saunders general form 
            C        O         V
        ----------
        |        |
     C  |   fc   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fc   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fc   |
        |        |        |        |
        ----------------------------
      */
        if (iter < scf_info->maxiter-1 && 
           !scf_info->converged && !scf_info->fock_typ) {
          if (occi == occj) 
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else if (occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if (occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
        }

     /* Guest & Saunders' form for high spin
            C        O         V
        ----------
        |        |
     C  | 2fc-fo |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo | 2fc-fo |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   | 2fc-fo |
        |        |        |        |
        ----------------------------
      */
        else if (iter < scf_info->maxiter-1 && !scf_info->converged 
                                          && scf_info->fock_typ == 1) {
          if ((occi == occj) || (occi && occj))
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if (occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
        }

     /* test form
            C        O         V
        ----------
        |        |
     C  |   fo   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fo   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fo   |
        |        |        |        |
        ----------------------------
      */
        else if (iter < scf_info->maxiter-1 && !scf_info->converged &&
                                              scf_info->fock_typ == 2) {
          if (occi == occj) 
            FOblk[i*jsz+j] = FOblk[i*jsz+j];
          else if (occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if (occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
        }

     /* form for converged wavefunction
            C        O         V
        ----------
        |        |
     C  |   fc   |
        |        |
        -------------------
        |        |        |
     O  | 2fc-fo |   fo   |
        |        |        |
        ----------------------------
        |        |        |        |
     V  |   fc   |   fo   |   fo   |
        |        |        |        |
        ----------------------------
      */
        else {
          if ((occi+occj)==4.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else if (occi == occj) 
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
          else if (occi && occj)
            Fblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if (occi==2.0 || occj==2.0)
            Fblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Fblk[i*jsz+j] = FOblk[i*jsz+j];
        }

        if (ib==jb && j==i) {
          if (occi == occ0 && occi) {
            Fblk[i*jsz+j] -= scf_info->lvl_shift;
          } else {
            if (occi) Fblk[i*jsz+j] -= 0.5*scf_info->lvl_shift;
          }
        }
      }
    }
  }
}

GLOBAL_FUNCTION double
scf_iter_elect_energy(scf_info,Fock,FockO)
scf_struct_t *scf_info;
dmt_matrix Fock;
dmt_matrix FockO;
{
  return
    scf_electronic_energy(scf_info,Hcore,Fock,Pmat,FockO,PmatO,SScr1,SScr2);
}

GLOBAL_FUNCTION double
scf_iter_delta_density(scf_info)
scf_struct_t *scf_info;
{
  return scf_rms_delta_density(scf_info,DPmat);
}

