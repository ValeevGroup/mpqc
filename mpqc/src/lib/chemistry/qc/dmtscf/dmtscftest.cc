
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <tmpl.h>

#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>

#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <math/dmt/matrix.h>

#include <chemistry/qc/intv2/int_libv2.h>
}

#include <util/keyval/ipv2c.h>

#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>

/*************************************************************************/
static void
init_mp()
{
  int nproc,me,host;
  int top,ord,dir;

  open0(&nproc,&me,&host);
  setarc0(&nproc,&top,&ord,&dir);
}

// this is done for libintv2 which still uses the old libipv2 stuff
static void 
init_ipv2()
{
  FILE *infile = fopen("mpqc.in","r");

  ip_initialize(infile,stdout);

  ip_cwk_add(":default");
  ip_cwk_add(":mpqc");

  ip_append_from_input("input",stdout);
  ip_cwk_clear();
  ip_cwk_add(":default");
  ip_cwk_add(":mpqc");
}

static void
init_dmt(centers_t *centers)
{
  int i;
  int nshell, nshtr;
  int *shellmap;

  int_initialize_1e(0,0,centers,centers);
  int_initialize_offsets1(centers,centers);

  nshell = centers->nshell;
  nshtr = nshell * (nshell+1) / 2;

  shellmap = (int *) malloc(sizeof(int)*nshell);
  bzero(shellmap,sizeof(int)*nshell);

  for (i=0; i < nshell; i++) shellmap[i] = INT_SH_NFUNC(centers,i);

  dmt_def_map2(centers->nfunc, centers->nshell, shellmap, 0, 1);

  free(shellmap);

  if (mynode0() == 0) printf("\n");
  dmt_map_examine();

  int_done_offsets1(centers,centers);
  int_done_1e();
}

int
main()
{
  int errcod;
  int save_vector=1;
  int save_fock=0;

  centers_t centers, oldcenters;
  scf_struct_t scf_info;
  sym_struct_t sym_info;

  dmt_matrix Scf_Vec, Fock, FockO;

  init_mp();

  if (mynode0() == 0) {
   // initialize the keyval ipv2 C functions
    init_ipv2();

   // initialize all of the structs needed by the SC libraries
    errcod = scf_init_scf(&centers, &scf_info, &sym_info, "mpqc.in");

   // pretty print the scf options
    scf_print_options(stdout, scf_info);

    ip_boolean("save_vector",&save_vector,0);
    ip_boolean("save_fock",&save_fock,0);
  }

  bcast0_scf_struct(&scf_info,0,0);
  bcast0_sym_struct(&sym_info,0,0);

  bcast0(&save_fock,sizeof(int),mtype_get(),0);
  bcast0(&save_vector,sizeof(int),mtype_get(),0);

  bcast0_centers(&centers,0,0);

 // if we need a projected guess, initialize oldcenters
  if (scf_info.proj_vector) {
    errcod = scf_make_old_centers(&centers, &oldcenters);
  }

 // free up the ipv2 memory
  ip_done();

 // initialize the dmt library
  init_dmt(&centers);

 // we need dmt matrices to hold the scf vector and Fock matrices
  Scf_Vec = dmt_create("scf vector",scf_info.nbfao,COLUMNS);
  Fock = dmt_create("fock matrix",scf_info.nbfao,SCATTERED);
  if (scf_info.iopen)
    FockO = dmt_create("open fock matrix",scf_info.nbfao,SCATTERED);
  else
    FockO = dmt_nil();

 // read the scf vector if it's there
  if (scf_info.restart) {
    char vecfile[512];
    sprintf(vecfile,"./%s.scfvec",scf_info.fname);
    dmt_read(vecfile,Scf_Vec);
  }

 // calculate the scf vector
  errcod = scf_vector(&scf_info, &sym_info, &centers, Fock, FockO, Scf_Vec,
                      &oldcenters, stdout);

 // save the scf vector
  if (save_vector) {
    char vecfile[512];
    sprintf(vecfile,"./%s.scfvec",scf_info.fname);
    dmt_write(vecfile,Scf_Vec);
  }

 // save the fock matrices
  if (save_fock) {
    char fockfile[512];
    sprintf(fockfile,"./%s.fock",scf_info.fname);
    dmt_write(fockfile,Fock);

    if (scf_info.iopen) {
      sprintf(fockfile,"./%s.fock",scf_info.fname);
      dmt_write(fockfile,Fock);
    }
  }

  exit(0);
}
