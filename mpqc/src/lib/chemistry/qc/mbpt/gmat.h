
#ifndef _chemistry_qc_mbpt_gmat_h
#define _chemistry_qc_mbpt_gmat_h
 
int
mbpt_make_gmat(scf_struct_t *scf_info, centers_t *centers,
          RefSymmSCMatrix& Gmat, double *DPmat, FILE *outfile);

int
mbpt_init_gmat(centers_t *centers, scf_struct_t *scf_info, double *intbuf);

void
mbpt_done_gmat(centers_t *centers, scf_struct_t *scf_info);

#endif
