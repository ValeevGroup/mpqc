
#include <stdio.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf.h>
         
#include <chemistry/qc/dmtscf/scf_mkden.gbl>
#include <chemistry/qc/dmtscf/scf_gmat.gbl>
#include <chemistry/qc/dmtscf/scf_oeis.gbl>

#include "scf_ffo.gbl"
#include "scf_ffo.lcl"

/* computes FOCK and FOCKO */
GLOBAL_FUNCTION int
scf_ffo(S, _scf_info, _sym_info, _centers, SCF_VEC, FOCK, FOCKO)
    dmt_matrix S;
    scf_struct_t *_scf_info;
    sym_struct_t *_sym_info;
    centers_t *_centers;
    dmt_matrix SCF_VEC;
    dmt_matrix FOCK;
    dmt_matrix FOCKO;
{
  double_vector_t occ_num;
  dmt_matrix GMAT;
  dmt_matrix GMATO;
  dmt_matrix PMAT;
  dmt_matrix DPMAT;
  dmt_matrix PMATO;
  dmt_matrix DPMATO;
  dmt_matrix SCR1;
  dmt_matrix SCR2;
  dmt_matrix SCR3;
  dmt_matrix SSCR1;
  dmt_matrix SSCR2;
  FILE *_outfile = stdout;
  int nbasis = _scf_info->nbfao;
  int errcod;
  int j;

/* set up the occupation number array */
  errcod = allocbn_double_vector(&occ_num,"n",nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"could not allocate memory for occ_num vector\n");
    return(-1);
    }

  for(j=0; j < _scf_info->nclosed ; j++)                occ_num.d[j]=2.0;
  for(; j < _scf_info->nclosed+_scf_info->nopen ; j++) 
                                                            occ_num.d[j]=1.0;
  for(; j < _scf_info->nbfso ; j++)                    occ_num.d[j]=0.0;


/* form initial density matrices */

  PMAT = dmt_create("libscfv3 density matrix",nbasis,SCATTERED);
  dmt_fill(PMAT,0.0);
  DPMAT = dmt_create("libscfv3 density diff matrix",nbasis,SCATTERED);
  dmt_fill(DPMAT,0.0);
  if(_scf_info->iopen) {
    PMATO = dmt_create("libscfv3 open density matrix",nbasis,SCATTERED);
    dmt_fill(PMATO,0.0);
    DPMATO = dmt_create("libscfv3 open density diff matrix",nbasis,SCATTERED);
    dmt_fill(DPMATO,0.0);
    }
  else {
    PMATO = dmt_nil();
    DPMATO = dmt_nil();
    }

/*  errcod = scf_make_density(_scf_info,_irreps,
                           SCF_VEC,PMAT,DPMAT,PMATO,DPMATO,&occ_num,_outfile);
    IMBN: replaced with the following two lines to conform with new 
          082294 version of scf */

  errcod = scf_make_density(_scf_info,
                           SCF_VEC,PMAT,DPMAT,PMATO,DPMATO,occ_num.d);

  if(errcod!=0) {
    fprintf(_outfile,"trouble forming density matrices\n");
    return(-1);
    }

/* and let's form the pk file */

  GMAT = dmt_create("libscfv3 g matrix",nbasis,SCATTERED);
  dmt_fill(GMAT,0.0);
  if(_scf_info->iopen) {
    GMATO = dmt_create("libscfv3 open g matrix",nbasis,SCATTERED);
    dmt_fill(GMATO,0.0);
    }
  else {
    GMATO = dmt_nil();
    }

/* allocate scratch arrays now */

  SCR1 = dmt_create("scf_iter: scr1",nbasis,COLUMNS);
  SCR2 = dmt_create("scf_iter: scr2",nbasis,COLUMNS);
  SCR3 = dmt_create("scf_iter: scr3",nbasis,COLUMNS);
  SSCR1 = dmt_create("scf_iter: scr4",nbasis,SCATTERED);
  SSCR2 = dmt_create("scf_iter: scr5",nbasis,SCATTERED);

/* form g matrix */
/*  errcod = scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,
                         GMAT,GMATO,DPMAT,DPMATO,SSCR1,SSCR2,1,0,_outfile);
    IMBN: replaced with the following two lines to conform with new
          082294 version of scf */

  scf_init_gmat(_centers,_scf_info);
  errcod = scf_make_gmat(_scf_info,_sym_info,_centers,
                         GMAT,GMATO,DPMAT,DPMATO,SSCR1,SSCR2,_outfile);
  scf_done_gmat(_centers,_scf_info);

  if(errcod != 0) {
    fprintf(_outfile,"scf_vector: trouble forming gmat\n");
    return(-1);
    }

  make_fock(S,_centers,_scf_info,_sym_info,FOCK,FOCKO,GMAT,GMATO,SSCR1,SSCR2,_outfile);

  dmt_free(PMAT);
  dmt_free(GMAT);
  dmt_free(DPMAT);
  dmt_free(SCR1);
  dmt_free(SCR2);
  dmt_free(SCR3);
  dmt_free(SSCR1);
  dmt_free(SSCR2);
  if(_scf_info->iopen) {
    dmt_free(PMATO);
    dmt_free(GMATO);
    dmt_free(DPMATO);
    }
  free_double_vector(&occ_num);

  return 0;
}

LOCAL_FUNCTION VOID
make_fock(S,_centers,_scf_info,_sym_info,FOCK,FOCKO,GMAT,GMATO,SSCR1,SSCR2,_outfile)
dmt_matrix S;
centers_t *_centers;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
dmt_matrix FOCK;
dmt_matrix FOCKO;
dmt_matrix GMAT;
dmt_matrix GMATO;
dmt_matrix SSCR1;
dmt_matrix SSCR2;
FILE *_outfile;
{
  int errcod;
  int iopen=_scf_info->iopen;
  int nbasis = _scf_info->nbfao;
  dmt_matrix T = dmt_create("libscfv3 kinetic matrix",nbasis,SCATTERED);
  dmt_matrix V = dmt_create("libscfv3 potential matrix",nbasis,SCATTERED);
  dmt_matrix HCORE = dmt_create("libscfv3 hcore matrix",nbasis,SCATTERED);

 /* form full g matrix from the skeleton gmatrix, place the result in SSCR1 */

  if(_sym_info->g > 1) {
    sym_sym_matrix(_centers,_sym_info,GMAT,SSCR1);
    if(iopen) sym_sym_matrix(_centers,_sym_info,GMATO,SSCR2);
    }
  else {
    dmt_copy(GMAT,SSCR1);
    if(iopen) dmt_copy(GMATO,SSCR2);
    }

  int_done_offsets2(_centers,_centers,_centers,_centers);
  int_done_erep();

  tim_enter("scf_oeis (again)");
/* errcod = scf_oeis(_scf_info,_sym_info,_irreps,_centers,S,T,V,HCORE,_outfile);
   IMBN: changed to following line to conform with new 082294 scf */
  errcod = scf_oeis(_scf_info,_centers,S,T,V,HCORE,_outfile);
  tim_exit("scf_oeis (again)");
  if(errcod != 0) {
    fprintf(stderr,"scf_ffo:\n");
    fprintf(stderr,"trouble in forming one-electron integrals\n");
    abort(); }
  dmt_free(T);
  dmt_free(V);

 /* F = H + G
  * FO = H + G - GO */
  dmt_copy(HCORE,FOCK);
  dmt_free(HCORE);
  dmt_sum(SSCR1,FOCK);

  if(iopen) {
    dmt_copy(SSCR2,FOCKO);
    dmt_scale(FOCKO,-1.0);
    dmt_sum(FOCK,FOCKO);
    }
  }
