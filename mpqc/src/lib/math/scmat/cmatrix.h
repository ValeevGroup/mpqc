
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
    void cmat_schmidt(double **rows, double *S, int nrow, int nc);


#ifdef __cplusplus
}
#endif

#endif
