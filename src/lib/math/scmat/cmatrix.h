/*
 * cmatrix.h
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

#ifndef _math_scmat_cmatrix_h
#define _math_scmat_cmatrix_h

#ifdef __cplusplus
extern "C" {
#endif
    
/* These routines work in terms of simple arrays.  No testing is done
 * on the input.
 */
    double cmat_determ(double**matrix,int sym,int dim);
    double cmat_invert(double**matrix,int sym,int dim);
    double cmat_solve_lin(double**,int sym,double*,int dim);
    void cmat_mxm(double**a,int transpose_a,
                  double**b,int transpose_b,
                  double**c,int transpose_c,
                  int nrow, int nlink, int ncol,
                  int add);
    void cmat_symmetric_mxm(double**a,int na, /* a is (na,na) */
                            double**b,int nb, /* b is (na,nb) */
                            int add);
    void cmat_transform_symmetric_matrix(double**a,int na, /* a is (na,na) */
                                         double**b,int nb, /* b is (nb,nb) */
                                         double**c,        /* c is (na,nb) */
                                         int add);
    void cmat_transform_diagonal_matrix(double**a,int na, /* a is (na,na) */
                                        double*b,int nb,  /* b is (nb,nb) */
                                        double**c,        /* c is (na,nb) */
                                        int add);
    double** cmat_new_square_matrix(int n);
    double** cmat_new_rect_matrix(int n,int m);
    void cmat_delete_matrix(double**matrix);
    void cmat_transpose_square_matrix(double**matrix,int n);
    void cmat_transpose_matrix(double**a,int nrow,int ncol);
    void cmat_matrix_pointers(double**ptrs,double*matrix,int nrow, int ncol);
    void cmat_diag(double**symm_a, double*evals, double**evecs, int n,
                   int matz, double tol);
    void cmat_eigensystem(/*const*/ double**symm_a, /*const*/ double**symm_s, double*evals, double**evecs, int n,
                          int matz);
    void cmat_schmidt(double **rows, double *S, int nrow, int nc);
    int cmat_schmidt_tol(double **C, double *S, int nrow, int ncol,
                         double tolerance, double *res);


#ifdef __cplusplus
}
#endif

#endif
