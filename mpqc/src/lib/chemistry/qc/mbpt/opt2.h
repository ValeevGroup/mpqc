
#ifndef _chemistry_qc_mbpt_opt2_h
#define _chemistry_qc_mbpt_opt2_h

int
mbpt_opt2(centers_t &centers, scf_struct_t &scf_info, sym_struct_t &sym_info,
          dmt_matrix Scf_Vec, dmt_matrix Fock, dmt_matrix FockO,
          int nfzc, int nfzv, int mem,
          int v1, int v2, int lb,
          FILE* outfile);

int
mbpt_opt2_v1(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
             double_vector_t *_evals, int nfzc, int nfzv, int mem,
             FILE* outfile);

int
mbpt_opt2_v2(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
             double_vector_t *_evals, int nfzc, int nfzv, int mem,
             FILE* outfile);

int
mbpt_opt2v2lb(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
              double_vector_t *_evals, int nfzc, int nfzv, int mem,
              FILE* outfile);

#endif
