
#ifndef _opt2_h
#define _opt2_h

int
opt2_v1(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
            double_vector_t *_evals, int nfzc, int nfzv, int mem,FILE* outfile);

int
opt2_v2(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
            double_vector_t *_evals, int nfzc, int nfzv, int mem,FILE* outfile);

int
opt2v2lb(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
            double_vector_t *_evals, int nfzc, int nfzv, int mem,FILE* outfile);

#endif
