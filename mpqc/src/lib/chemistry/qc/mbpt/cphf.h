#ifndef _mpqcic_cphf_h
#define _mpqcic_cphf_h
 
void
mbpt_cphf(centers_t *centers, scf_struct_t *scf_info, FILE *outfile,
          int nbasis, int nvir, int nocc, double **Scf_Vec, double *L,
          double *eigval, RefSCMatrix& P2);

#endif
