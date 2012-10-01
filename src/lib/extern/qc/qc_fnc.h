//
// qc_fnc.h
//
// Copyright (C) 2012 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUG__
#pragma interface
#endif

// wraps most interesting MPQC classes for use via simple function/subroutine calls
// since this is meant to be used from low-level languages (C, F*) no namespaces, and C-linkage

#ifndef _mpqc_src_lib_extern_qc_qcfnc_h
#define _mpqc_src_lib_extern_qc_qcfnc_h

#if __cplusplus
extern "C" {
#endif

/**
 * initialize a molecule data using atomic numbers and Cartesian coordinates
 * @param[in] natoms the number of atoms
 * @param[in] Z atomic numbers of the atoms
 * @param[in] xyz Cartesian coordinates of the atoms in atomic units (row-major matrix with natoms rows and 3 columns)
 * @param[in] use_symmetry if true, will detect the point group automatically
 */
void init_molecule_(int natoms, const double* Z, const double* xyz, int use_symmetry);
/**
 * initialize a molecule using a file in XYZ format
 * @param[in] fname file name
 * @param[in] use_symmetry if true, will detect the point group automatically
 * @param[in] fname_nchar number of characters in fname, this parameter is implicit when calling from Fortran, hence no need to provide it
 */
void init_molecule_xyz_(const char* fname, int use_symmetry, int fname_nchar);

/**
 * initialize the orbital basis set
 * @param[in] basis_name the basis set name
 * @param[in] basis_name_nchar number of characters in basis_name, this parameter is implicit when calling from Fortran, hence no need to provide it
 */
void init_basis_set_(const char* basis_name, int basis_name_nchar);

/**
 * initialize the orbital basis set using a file in Gaussian 94 format
 * @param[in] fname the file name with the basis set specification
 * @param[in] fname_nchar number of characters in fname, this parameter is implicit when calling from Fortran, hence no need to provide it
 */
void init_basis_set_g94_(const char* fname, int fname_nchar);

/**
 * @return the basis set size
 */
int basis_set_nbasis_();

/**
 * @return the number of shells in the basis
 */
int basis_set_nshell_();

/**
 * @return the maximum number of basis functions in a shell of this basis
 */
int basis_set_max_nfunction_in_shell_();

/**
 * @return the array that contains the number of basis functions in each shell
 */
int* basis_set_shell_to_nfunction_();

/**
 * subroutine version of basis_set_shell_to_nfunction_(), to be used from pre-2003 fortran
 * @param[out] s2nf the array that contains the number of basis functions in each shell
 */
void basis_set_shell_to_nfunction_subrt_(int* s2nf);

/**
 * @return the array that contains the index of the first basis function from each shell
 */
int* basis_set_shell_to_function_();

/**
 * subroutine version of basis_set_shell_to_function_(), to be used from pre-2003 fortran
 * @param[out] s2f the array that contains the index of the first basis function from each shell
 */
void basis_set_shell_to_function_subrt_(int* s2f);

/**
 * @return the array that contains the center (atom) on which each shell is centered
 */
int* basis_set_shell_to_center_();

/**
 * subroutine version of basis_set_shell_to_center_(), to be used from pre-2003 fortran
 * @param[out] s2c the array that contains the center (atom) on which each shell is centered
 */
void basis_set_shell_to_center_subrt_(int* s2c);

/**
 * initialize the integrals factory
 */
void init_integrals_();
/**
 * de-initialize the integrals factory
 */
void done_integrals_();

/**
 * initialize the overlap integral evaluator
 */
void init_overlap_integrals_();
/**
 * if you are using Fortran older than 2003 use set_overlap_integrals_buffer_() instead
 * @return the pointer to the buffer that will contain overlap integrals
 */
const double* overlap_integrals_buffer_();
/**
 * use this instead of overlap_integrals_buffer_() to provide an external buffer to store integrals in
 * @param[in] buf the pointer to the buffer that will contain overlap integrals
 */
void set_overlap_integrals_buffer_(double* buf);
/**
 * compute a shell-pair of overlap integrals
 * @param[in] bra bra shell index
 * @param[in] ket ket shell index
 */
void compute_overlap_shell_(int bra, int ket);
/**
 * de-initialize the overlap integral evaluator
 */
void done_overlap_integrals_();

/**
 * initialize the core Hamiltonian (T+V1) integral evaluator
 */
void init_hcore_integrals_();
/**
 * if you are using Fortran older than 2003 use set_hcore_integrals_buffer_() instead
 * @return the pointer to the buffer that will contain core Hamiltonian integrals
 */
const double* hcore_integrals_buffer_();
/**
 * use this instead of hcore_integrals_buffer_() to provide an external buffer to store integrals in
 * @param[in] buf the pointer to the buffer that will contain core Hamiltonian integrals
 */
void set_hcore_integrals_buffer_(double* buf);
/**
 * compute a shell-pair of core Hamiltonian integrals
 * @param[in] bra bra shell index
 * @param[in] ket ket shell index
 */
void compute_hcore_shell_(int bra, int ket);
/**
 * de-initialize the core Hamiltonian integral evaluator
 */
void done_hcore_integrals_();

/**
 * initialize the 2-e Coulomb (electron repulsion) integral evaluator
 */
void init_twoecoulomb_integrals_();
/**
 * if you are using Fortran older than 2003 use set_twoecoulomb_integrals_buffer_() instead
 * @return the pointer to the buffer that will contain 2-e Coulomb integrals
 */
const double* twoecoulomb_integrals_buffer_();
/**
 * use this instead of twoecoulomb_integrals_buffer_() to provide an external buffer to store integrals in
 * @param[in] buf the pointer to the buffer that will contain 2-e Coulomb integrals
 */
void set_twoecoulomb_integrals_buffer_(double* buf);
/**
 * compute a shell-pair of 2-e Coulomb integrals, (bra1 ket1|bra2 ket2)
 * @param[in] bra1 bra1 shell index
 * @param[in] ket1 ket1 shell index
 * @param[in] bra2 bra2 shell index
 * @param[in] ket2 ket2 shell index
 */
void compute_twoecoulomb_shell_(int bra1, int ket1, int bra2, int ket2);
/**
 * de-initialize the 2-e Coulomb integral evaluator
 */
void done_twoecoulomb_integrals_();

#if __cplusplus
}; // end of extern "C"
#endif

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
