
/* This computes the exchange energy of the wavefunction. */

/* $Log$
 * Revision 1.2  1993/12/30 13:31:17  etseidl
 * merge in clj changes, do global sum of exchange energy in scf_ex.c
 *
 * Revision 1.1.1.1  1993/12/29  12:53:15  etseidl
 * SC source tree 0.1
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <math/array/math_lib.h>
#include <util/misc/libmisc.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#include "scf_ex.gbl"
#include "scf_ex.lcl"

GLOBAL_FUNCTION void
scf_ex(_scf_info,_centers,PMAT)
scf_struct_t *_scf_info;
centers_t*_centers;
dmt_matrix PMAT;
{
  loop_t*loop;
  int ish,jsh,ksh,lsh, isz,jsz,ksz,lsz;
  int i,j,k,l;
  int nl;
  double *pkl, *pij;
  double *buf;
  double exchange1 = 0.0;
  double exchange2 = 0.0;
  int ij,kl;
  int index;
  int nlocal=dmt_nlocal(PMAT);
  int flags;

/*
  dmt_printf(" %6.4f",PMAT);
 */

  if (_scf_info->iopen) {
      fprintf(stderr,"scf_ex: can't compute exchange energy for open shells\n");
      return;
    }

  tim_enter("scf_ex");

  /* initialize the integrals routines */
  int_initialize_offsets2(_centers,_centers,_centers,_centers);
  flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  buf = int_initialize_erep(flags,0,_centers,_centers,_centers,_centers);

  loop = dmt_ngl_create("%m",PMAT);
  while(dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&ksh,&ksz,&lsh,&lsz,&pkl)) {
      double factor,factor2;

      for(nl=0; nl < nlocal; nl++) {
        dmt_get_block_dsc(PMAT,nl,&ish,&isz,&jsh,&jsz,&pij);

	if (ish != jsh && ksh != lsh) factor = 2.0;
	else factor = 1.0;

	if (ish != jsh) factor *= 0.5;
        if (ksh != lsh) factor *= 0.5;

	/* compute the integrals for this shell block */
	int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&ish,&ksh,&jsh,&lsh);
	/* printf("(%d %d|%d %d)[0] is %12.8f\n",ish,ksh,jsh,lsh,buf[0]); */

        index = 0;
	for (i=0; i<isz; i++) {
	    for (k=0; k<ksz; k++) {
		for (j=0; j<jsz; j++) {
                    ij = i*jsz + j;
		    for (l=0; l<lsz; l++) {
                        kl = k*lsz + l;

			factor2 = factor;
			if (ish == jsh && i != j) factor2 *= 0.5;
			if (ksh == lsh && k != l) factor2 *= 0.5;

			exchange1 += buf[index]*pij[ij]*pkl[kl]*factor2;
/*
			printf("adding(1) P%d%d P%d%d (%d%d|%d%d) % 6.4f * % 6.4f * % 6.4f * %2.1f\n",
			       ish,jsh,ksh,lsh,ish,ksh,jsh,lsh,
			       buf[index],pij[ij],pkl[kl],factor);
 */
			index++;
		      }
		  }
	      }
	  }

        if (ish == jsh && ksh == lsh) continue;

	int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&ish,&lsh,&ksh,&jsh);
	/* printf("(%d %d|%d %d)[0] is %12.8f\n",ish,lsh,ksh,jsh,buf[0]); */

        index = 0;
	for (i=0; i<isz; i++) {
	    for (l=0; l<lsz; l++) {
		for (k=0; k<ksz; k++) {
                    kl = k*lsz + l;
		    for (j=0; j<jsz; j++) {
                        ij = i*jsz + j;

			factor2 = factor;
			if (ish == jsh && i != j) factor2 *= 0.5;
			if (ksh == lsh && k != l) factor2 *= 0.5;

			exchange2 += buf[index]*pij[ij]*pkl[kl]*factor2;
/*
			printf("adding(2) P%d%d P%d%d (%d%d|%d%d) % 6.4f * % 6.4f * % 6.4f * %2.1f\n",
			       ish,jsh,ksh,lsh,ish,lsh,ksh,jsh,
			       buf[index],pij[ij],pkl[kl],factor);
 */
			index++;
		      }
		  }
	      }
	  }
      }

    }
  }

/* we are finished with loop */
  dmt_ngl_kill(loop);

  int_done_erep();
  int_done_offsets2(_centers,_centers,_centers,_centers);

/*
  printf("exchange1 is %14.8f\n",-0.25*exchange1);
  printf("exchange2 is %14.8f\n",-0.25*exchange2);
 */

  gsum0(&exchange1,1,5,mtype_get(),0);
  gsum0(&exchange2,1,5,mtype_get(),0);

  if (mynode0()==0)
    printf("\n  The exchange energy is %14.8f\n",-0.25*(exchange1+exchange2));

  tim_exit("scf_ex");
}
