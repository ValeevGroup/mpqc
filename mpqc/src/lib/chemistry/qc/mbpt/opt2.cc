
#include <stdio.h>
#include <math.h>

#include <util/group/picl.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/mbpt/ffo.h>
#include <chemistry/qc/mbpt/opt2.h>

#include <chemistry/qc/mbpt/opt2_fock.h>

int
mbpt_opt2(centers_t &centers, scf_struct_t &scf_info, sym_struct_t &sym_info,
          dmt_matrix Scf_Vec, dmt_matrix Fock, dmt_matrix FockO,
          int nfzc, int nfzv, int mem_alloc,
          int do_opt2_v1, int do_opt2_v2, int do_opt2v2lb,
          FILE* outfile)
{
  dmt_matrix S = dmt_create("libscfv3 overlap matrix",
                            scf_info.nbfao,SCATTERED);
  dmt_matrix SAHALF;
  dmt_matrix SC;
  dmt_matrix EVECS;
  dmt_matrix SCR;
  double_vector_t occ_num;
  double_vector_t evals;
  dmt_matrix SCR1,SCR2,SCR3;
  // this got free'ed somewhere
  if(scf_info.iopen)
      FockO = dmt_create("opt2 open fock matrix",scf_info.nbfao,SCATTERED);
  SCR1 = dmt_create("opt2:scr1",scf_info.nbfao,COLUMNS);
  SCR2 = dmt_create("opt2:scr2",scf_info.nbfao,COLUMNS);
  SCR3 = dmt_create("opt2:scr3",scf_info.nbfao,COLUMNS);

  if (!do_opt2_v1 && !do_opt2_v2) do_opt2_v2 = 1;

  tim_enter("opt2");

  mbpt_ffo(S, &scf_info, &sym_info, &centers, Scf_Vec, Fock, FockO);

  if (scf_info.iopen) {
      mbpt_make_opt2_fock(&scf_info,Fock,FockO,Scf_Vec,SCR1,SCR2,SCR3);
      dmt_free(SCR3);
      dmt_copy(Fock,SCR1); /* dmt_diag needs a columns dist. matrix */
      allocbn_double_vector(&evals,"n",scf_info.nbfao);
      dmt_diag(SCR1,SCR2,evals.d); /*SCR2 transforms from old to new mo basis */
      dmt_copy(Scf_Vec,SCR1);
      dmt_transpose(SCR1);
      dmt_mult(SCR1,SCR2,Scf_Vec);
      dmt_free(SCR1);
      dmt_free(SCR2);
    }
  else {
      /* form 'sahalf' matrix sahalf = u*ei^-0.5*u~ */
      allocbn_double_vector(&evals,"n",scf_info.nbfao);
      SAHALF= dmt_create("libscfv3 scf_core_guess scr4",
                         scf_info.nbfao,COLUMNS);
      EVECS = dmt_create("libscfv3 scf_core_guess scr3",
                         scf_info.nbfao,COLUMNS);
      SCR = dmt_create("libscfv3 scf_core_guess scr5",scf_info.nbfao,COLUMNS);
      SC = dmt_columns("libscfv3 scf_core_guess scr1",S);
      /* diagonalize overlap matrix */
      dmt_diag(SC,EVECS,evals.d);
      /* form SAHALF matrix (s^(-1/2), Sz&Ostl p. 143) */
      for(int i=0; i < scf_info.nbfao; i++) evals.d[i] = 1.0/sqrt(evals.d[i]);
      dmt_fill(SAHALF,0.0);
      dmt_set_diagonal(SAHALF,evals.d);
      /* form the orthogonalization matrix S^(-1/2) (Szabo&Ostlund p. 143)
       * (called SAHALF here) */
      dmt_transpose(EVECS);
      dmt_mult(SAHALF,EVECS,SCR);
      dmt_mult(EVECS,SCR,SAHALF);
      dmt_free(EVECS);
      dmt_free(SC);
      dmt_free(SCR);
      dmt_free(S);
      dmt_copy(Fock,SCR1); /* need a columns distr. matrix */
      dmt_mult(SCR1,SAHALF,SCR2);
      dmt_mult(SAHALF,SCR2,SCR3);
      /* SCR3 is now the Fock matrix in the orthogonalized ao basis */
      dmt_diag(SCR3,Scf_Vec,evals.d);
      dmt_copy(Scf_Vec,SCR1);
      dmt_mult(SAHALF,SCR1,Scf_Vec); /* Sz&Ostl p.146 point 9 */
      dmt_free(SCR1);
      dmt_free(SCR2);
      dmt_free(SCR3);
      }

    if (do_opt2_v1) {
        sync0();
        tim_enter("opt2_v1");
        mbpt_opt2_v1(&centers,&scf_info,Scf_Vec,&evals,nfzc,nfzv,mem_alloc,
                     outfile);
        tim_exit("opt2_v1");
      }
    if (do_opt2_v2) {
        sync0();
        tim_enter("opt2_v2");
        mbpt_opt2_v2(&centers,&scf_info,Scf_Vec,&evals,nfzc,nfzv,mem_alloc,
                     outfile);
        tim_exit("opt2_v2");
      }
    if (do_opt2v2lb) {
        sync0();
        tim_enter("opt2v2lb");
        mbpt_opt2v2lb(&centers,&scf_info,Scf_Vec,&evals,nfzc,nfzv,mem_alloc,
                      outfile);
        tim_exit("opt2v2lb");
      }

    free_double_vector(&evals);

    tim_exit("opt2");
    }
