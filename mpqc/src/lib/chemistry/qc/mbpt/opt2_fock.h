
#ifndef _chemistry_qc_mbpt_opt2_fock_h
#define _chemistry_qc_mbpt_opt2_fock_h

#ifdef __cplusplus
extern "C" {
#endif

int
mbpt_make_opt2_fock(scf_struct_t *scf_info,
                    dmt_matrix FOCK, dmt_matrix FOCKO,
                    dmt_matrix SCF_VEC,
                    dmt_matrix SCR1, dmt_matrix SCR2, dmt_matrix SCR3);

#ifdef __cplusplus
}
#endif

#endif
