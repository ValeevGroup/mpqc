
#ifndef _chemistry_qc_intv2_h
#define _chemistry_qc_intv2_h

#include <chemistry/qc/intv2/atoms.h>

/* integrals and target may overlap */
void int_transform_1e(double *integrals, double *target,
                      shell_t *sh1, shell_t *sh2);

/* integrals and target may not overlap */
void int_accum_transform_1e(double *integrals, double *target,
                            shell_t *sh1, shell_t *sh2);

/* integrals and target may overlap */
void int_transform_1e_xyz(double *integrals, double *target,
                          shell_t *sh1, shell_t *sh2);

/* integrals and target may not overlap */
void int_accum_transform_1e_xyz(double *integrals, double *target,
                                shell_t *sh1, shell_t *sh2);

/* integrals and target may overlap */
void int_transform_2e(double *integrals, double *target,
                      shell_t *sh1, shell_t *sh2,
                      shell_t *sh3, shell_t *sh4);

#endif
