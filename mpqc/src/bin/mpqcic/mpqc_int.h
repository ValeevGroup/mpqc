

/* mpqc_int.h -- header for mpqc interface routines
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      March, 1993
 */

#ifndef _libGeom_mpqc_int_h
#define _libGeom_mpqc_int_h

enum geom_codes { GEOM_DONE, GEOM_ABORT, GEOM_COMPUTE_ENERGY,
                  GEOM_COMPUTE_GRADIENT, GEOM_COMPUTE_HESSIAN };

#include <stdio.h>

extern "C" {
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math/dmt/libdmt.h>
}

#include <util/keyval/keyval.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>


RefSCDimension Geom_dim_natom3();

int Geom_init_mpqc(RefMolecule&, const RefKeyVal&, const char *);

int Geom_update_mpqc(double, RefSCVector&, const RefKeyVal&, const char *);

void Geom_done_mpqc(const RefKeyVal&, int converged, const char *);
void Geom_write_pdb(const RefKeyVal&, RefMolecule&, const char *, char * =0);

int Scf_charges_from_esp(centers_t*,int,double_vector_t*,
                         double_vector_t*, expts_t*, double, int, FILE*,
                         const RefKeyVal&);

int do_mp2(centers_t&, scf_struct_t&, dmt_matrix, dmt_matrix, FILE*,
            const RefKeyVal&);


#endif
