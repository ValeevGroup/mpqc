#ifdef __cplusplus
extern "C" {
#endif

#ifndef __CBLAS_H__
#define __CBLAS_H__

/*   C BLAS  TYPE DEFINITIONS  */

typedef int  Integer;

typedef struct { 
            float real;
            float imag; 
} Complex;

typedef struct {
	    double real;
            double imag; 
} Zomplex;

typedef enum { NoTranspose,
               Transpose,
               ConjugateTranspose } MatrixTranspose;

typedef enum { UpperTriangle,
               LowerTriangle } MatrixTriangle;

typedef enum { UnitTriangular,
               NotUnitTriangular } MatrixUnitTriangular;

typedef enum { LeftSide,
               RightSide } OperationSide;


/********                                                             ********
 ********                    LEVEL 1 BLAS                             ********
 ********                                                             ********/

/*                         Generate a plane rotation                         */

extern void srotg(  float *a,  float *b,  float *c,  float *s );
extern void drotg( double *a, double *b, double *c, double *s );


/*                    Generate a modified plane rotation                     */

/*
extern void srotmg(  float *d1,  float *d2,  float *a,  float b, 
                     float *param );
extern void drotmg( double *d1, double *d2, double *a, double b, 
                    double *param );
*/


/*                          Apply a plane rotation                           */

extern void srot( Integer n,  float *x, Integer incx,  float *y,
	          Integer incy,  float c,  float s );
extern void drot( Integer n, double *x, Integer incx, double *y,
	          Integer incy, double c, double s );


/*                      Apply a modified plane rotation                      */

/*
extern void srotm( Integer n,  float *x, Integer incx,  float *y,
	           Integer incy,  float *param );
extern void drotm( Integer n, double *x, Integer incx, double *y,
	           Integer incy, double *param );
*/


/*                             Swap two vectors                              *
 *                                x <-> y                                    */

extern void sswap( Integer n,   float *x, Integer incx,   float *y, 
                   Integer incy );
extern void dswap( Integer n,  double *x, Integer incx,  double *y, 
                   Integer incy );
extern void cswap( Integer n, Complex *x, Integer incx, Complex *y, 
                   Integer incy );
extern void zswap( Integer n, Zomplex *x, Integer incx, Zomplex *y, 
                   Integer incy );


/*                              Scale a vector                               *
 *                               x <- alpha*x                                */

extern void  sscal( Integer n,   float alpha,   float *x, Integer incx );
extern void  dscal( Integer n,  double alpha,  double *x, Integer incx );
extern void  cscal( Integer n, Complex alpha, Complex *x, Integer incx );
extern void  zscal( Integer n, Zomplex alpha, Zomplex *x, Integer incx );
extern void csscal( Integer n,   float alpha, Complex *x, Integer incx );
extern void zdscal( Integer n,  double alpha, Zomplex *x, Integer incx );


/*                         Copy one vector to another                        *
 *                                  y <- x                                   */

extern void scopy( Integer n,   float *x, Integer incx,   float *y, 
                   Integer incy );
extern void dcopy( Integer n,  double *x, Integer incx,  double *y, 
                   Integer incy );
extern void ccopy( Integer n, Complex *x, Integer incx, Complex *y, 
                   Integer incy );
extern void zcopy( Integer n, Zomplex *x, Integer incx, Zomplex *y, 
                   Integer incy );


/*                 Scale a vector then add to another vector                 *
 *                             y <- alpha*x + y                              */

extern void saxpy( Integer n,   float alpha,   float *x, Integer incx, 
                     float *y, Integer incy );
extern void daxpy( Integer n,  double alpha,  double *x, Integer incx, 
                    double *y, Integer incy );
extern void caxpy( Integer n, Complex alpha, Complex *x, Integer incx, 
                   Complex *y, Integer incy );
extern void zaxpy( Integer n, Zomplex alpha, Zomplex *x, Integer incx, 
                   Zomplex *y, Integer incy );


/*                         Dot product of two vectors                        * 
 *                                 dot <- xTy                                */

extern float   sdot( Integer n,  float *x, Integer incx,
                                 float *y, Integer incy );
/*
extern double dsdot( Integer n,  float *x, Integer incx,
                                 float *y, Integer incy );
*/
extern double  ddot( Integer n, double *x, Integer incx,
                                double *y, Integer incy );


/*                  Dot product of two vectors plus a scalar                 * 
 *                             dot <- alpha + xTy                            */
/*
float sdsdot( Integer n, float alpha, float *x, Integer incx, float *y,
	      Integer incy );
*/


/*                         Dot product of two vectors                        * 
 *                                dotu <- xTy                                */

extern Complex cdotu( Integer n, Complex *x, Integer incx, 
                                 Complex *y, Integer incy );
extern Zomplex zdotu( Integer n, Zomplex *x, Integer incx, 
                                 Zomplex *y, Integer incy );


/*                    Dot conjugate product of two vectors                   * 
 *                                dotc <- xHy                                */

extern Complex cdotc( Integer n, Complex *x, Integer incx, 
                                 Complex *y, Integer incy );
extern Zomplex zdotc( Integer n, Zomplex *x, Integer incx, 
                                 Zomplex *y, Integer incy );


/*                            2-Norm of a vector                             *
 *                              nrm2 <- ||x||2                               */

extern float   snrm2( Integer n,   float *x, Integer incx );
extern double  dnrm2( Integer n,  double *x, Integer incx );
extern float  scnrm2( Integer n, Complex *x, Integer incx );
extern double dznrm2( Integer n, Zomplex *x, Integer incx );


/*                            1-Norm of a vector                             *
 *                              asum <- ||x||1                               */

extern float  sasum( Integer n,   float *x, Integer incx );
extern double dasum( Integer n,  double *x, Integer incx );


/*                         1-Norm of a complex vector                        *
 *                      asum <- ||Re(x)||1 + ||Im(x)||1                      */

extern float  scasum( Integer n, Complex *x, Integer incx );
extern double zdasum( Integer n, Zomplex *x, Integer incx );


/*            Index that will give the infinity-Norm of a vector             *
 *               amax <- first k such that |x(k)| = max|x(i)|                */

/*extern Integer isamax( Integer n,   float *x, Integer incx );*/
/*extern Integer idamax( Integer n,  double *x, Integer incx );*/

/*     amax <- first k such that                                             *
 *     |Re( x(k) )| + |Im( x(k) )| = max( |Re( x(k) )| + |Im( x(k) )| )      */

extern Integer icamax( Integer n, Complex *x, Integer incx );
extern Integer izamax( Integer n, Zomplex *x, Integer incx );



/********                                                             ********
 ********                    LEVEL 2 BLAS                             ********
 ********                                                             ********/

/*                      Matrix - vector multiplication                       *
 *                                                                           *
 *                        y <- alpha*A*x    + beta*y                         *
 *                        y <- alpha*A**T*x + beta*y           A is m-by-n   *
 *                        y <- alpha*A**H*x + beta*y                         */

extern void sgemv( MatrixTranspose trans, Integer m, Integer n,
	             float alpha,   float *a, Integer tda, 
                     float *x,    Integer incx,
                     float  beta,   float *y, Integer incy );
extern void dgemv( MatrixTranspose trans, Integer m, Integer n,
	            double alpha,  double *a, Integer tda, 
                    double *x,    Integer incx, 
                    double beta,   double *y, Integer incy );
extern void cgemv( MatrixTranspose trans, Integer m, Integer n,
	           Complex alpha, Complex *a, Integer tda, 
                   Complex *x,    Integer incx, 
                   Complex beta,  Complex *y, Integer incy );
extern void zgemv( MatrixTranspose trans, Integer m, Integer n,
	           Zomplex alpha, Zomplex *a, Integer tda, 
                   Zomplex *x,    Integer incx, 
                   Zomplex beta,  Zomplex *y, Integer incy );


/*                   Banded Matrix - vector multiplication                   *
 *                                                                           *
 *                        y <- alpha*A*x    + beta*y                         *
 *                        y <- alpha*A**T*x + beta*y       A is m-by-n band  *
 *                        y <- alpha*A**H*x + beta*y                         */

extern void sgbmv( MatrixTranspose trans,  
                   Integer m,     Integer n,  Integer kl, Integer ku,
                     float alpha,   float *a, Integer tda,
	             float *x,    Integer incx,
                     float beta,    float *y, Integer incy );
extern void dgbmv( MatrixTranspose trans,  
                   Integer m,     Integer n,  Integer kl, Integer ku,
                    double alpha,  double *a, Integer tda,
	            double *x,    Integer incx,
                    double beta,   double *y, Integer incy );
extern void cgbmv( MatrixTranspose trans,
                   Integer m,     Integer n,  Integer kl, Integer ku,
                   Complex alpha, Complex *a, Integer tda,
	           Complex *x,    Integer incx, 
                   Complex beta,  Complex *y, Integer incy );
extern void zgbmv( MatrixTranspose trans, 
                   Integer m,     Integer n,  Integer kl, Integer ku,
                   Zomplex alpha, Zomplex *a, Integer tda,
	           Zomplex *x,    Integer incx, 
                   Zomplex beta,  Zomplex *y, Integer incy );


/*                  Hermitian Matrix - vector multiplication                 *
 *                                                                           *
 *                        y <- alpha*A*x + beta*y      A is n-by-n Hermitian */

extern void chemv( MatrixTriangle uplo, Integer n,   Complex alpha,
	                  Complex *a,   Integer tda, Complex *x, Integer incx,
	                  Complex beta, Complex *y,  Integer incy );
extern void zhemv( MatrixTriangle uplo, Integer n,   Zomplex alpha,
	                  Zomplex *a,   Integer tda, Zomplex *x, Integer incx,
	                  Zomplex beta, Zomplex *y,  Integer incy );


/*              Hermitian Banded Matrix - vector multiplication              *
 *                                                                           *
 *                        y <- alpha*A*x + beta*y      A is n-by-n Hermitian *
 *                                                                 band      */

extern void chbmv( MatrixTriangle uplo, Integer n,   Integer k, Complex alpha,
	                  Complex *a,   Integer tda, Complex *x, Integer incx,
	                  Complex beta, Complex *y,  Integer incy );
extern void zhbmv( MatrixTriangle uplo, Integer n,   Integer k, Zomplex alpha,
	                  Zomplex *a,   Integer tda, Zomplex *x, Integer incx,
	                  Zomplex beta, Zomplex *y,  Integer incy );


/*              Hermitian Packed Matrix - vector multiplication              *
 *                                                                           *
 *                        y <- alpha*A*x + beta*y      A is n-by-n Hermitian *
 *                                                            in packed form */

extern void chpmv( MatrixTriangle uplo, Integer n, Complex alpha,
	                  Complex *ap,  Complex *x, Integer incx,
	                  Complex beta, Complex *y, Integer incy );
extern void zhpmv( MatrixTriangle uplo, Integer n, Zomplex alpha,
	                  Zomplex *ap,  Zomplex *x, Integer incx,
	                  Zomplex beta, Zomplex *y, Integer incy );


/*                  Symmetric Matrix - vector multiplication                 *
 *                                                                           *
 *                        y <- alpha*A*x + beta*y      A is n-by-n Symmetric */

extern void ssymv( MatrixTriangle uplo, 
		   Integer n,  float alpha,  float *a, Integer tda,
		    float *x, Integer incx,  float beta,
		    float *y, Integer incy );
extern void dsymv( MatrixTriangle uplo, 
		   Integer n, double alpha, double *a, Integer tda, 
		   double *x, Integer incx, double beta,
	           double *y, Integer incy );


/*              Symmetric Banded Matrix - vector multiplication              *
 *                                                                           *
 *                        y <- alpha*A*x + beta*y      A is n-by-n Symmetric *
 *                                                                 Band      */

extern void ssbmv( MatrixTriangle uplo, Integer n, Integer k,
		    float alpha,   float *a,  Integer tda,\
		    float *x,    Integer incx,  float beta,
		    float *y,    Integer incy );
extern void dsbmv( MatrixTriangle uplo, Integer n, Integer k,
		   double alpha,  double *a,  Integer tda, 
		   double *x,    Integer incx, double beta,
		   double *y,    Integer incy );


/*              Symmetric Packed Matrix - vector multiplication              *
 *                                                                           *
 *                        y <- alpha*A*x + beta*y      A is n-by-n Symmetric *
 *                                                                 Packed    */

extern void sspmv( MatrixTriangle uplo, Integer n,
		    float alpha,  float *ap,    float *x, Integer incx,
	             float beta,  float *y, Integer incy );
extern void dspmv( MatrixTriangle uplo, Integer n,
		   double alpha, double *ap,   double *x, Integer incx,
	            double beta, double *y, Integer incy );


/*                 Triangular Matrix - vector multiplication                 *
 *                                                                           *
 *                                y <- A*x                                   *
 *                                y <- A**T*x        A is n-by-n triangular  *
 *                                y <- A**H*x                                */

extern void strmv( MatrixTriangle uplo,       MatrixTranspose trans, 
	           MatrixUnitTriangular diag, Integer n,
		     float *a, Integer tda,    float *x, Integer incx );
extern void dtrmv( MatrixTriangle uplo, MatrixTranspose trans, 
	           MatrixUnitTriangular diag, Integer n,
		    double *a, Integer tda,   double *x, Integer incx );
extern void ctrmv( MatrixTriangle uplo, MatrixTranspose trans, 
	           MatrixUnitTriangular diag, Integer n, 
		   Complex *a, Integer tda,  Complex *x, Integer incx );
extern void ztrmv( MatrixTriangle uplo,       MatrixTranspose trans, 
	           MatrixUnitTriangular diag, Integer n, 
		   Zomplex *a, Integer tda,  Zomplex *x, Integer incx );


/*              Triangular Band Matrix - vector multiplication               *
 *                                                                           *
 *                                y <- A*x                                   *
 *                                y <- A**T*x        A is n-by-n triangular  *
 *                                y <- A**H*x                    band        */

extern void stbmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n, Integer k,
	      float *a, Integer tda,   float *x, Integer incx );
extern void dtbmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n, Integer k,
	     double *a, Integer tda,  double *x, Integer incx );
extern void ctbmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n, Integer k,
	    Complex *a, Integer tda, Complex *x, Integer incx );
extern void ztbmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n, Integer k,
	    Zomplex *a, Integer tda, Zomplex *x, Integer incx );


/*             Triangular Packed Matrix - vector multiplication              *
 *                                                                           *
 *                                y <- A*x                                   *
 *                                y <- A**T*x        A is n-by-n triangular  *
 *                                y <- A**H*x                    packed      */

extern void stpmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n,
	      float *ap,    float *x,  Integer incx );
extern void dtpmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n,
	     double *ap,   double *x,  Integer incx );
extern void ctpmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n,
	    Complex *ap,  Complex *x,  Integer incx );
extern void ztpmv( MatrixTriangle uplo, MatrixTranspose trans, 
	    MatrixUnitTriangular diag, Integer n,
	    Zomplex *ap,  Zomplex *x,  Integer incx );


/*                  Solve a triangular system of equations                   *
 *                                                                           *
 *                             y <- A**(-1)*x                                *
 *                             y <- A**(-T)*x        A is n-by-n triangular  *
 *                             y <- A**(-H)*x                                */

extern void strsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	      float *a, Integer tda,   float *x, Integer incx );
extern void dtrsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	     double *a, Integer tda,  double *x, Integer incx );
extern void ctrsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	    Complex *a, Integer tda, Complex *x, Integer incx );
extern void ztrsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	    Zomplex *a, Integer tda, Zomplex *x, Integer incx );


/*               Solve a banded triangular system of equations               *
 *                                                                           *
 *                             y <- A**(-1)*x                                *
 *                             y <- A**(-T)*x        A is n-by-n triangular  *
 *                             y <- A**(-H)*x                    band        */

extern void stbsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n, Integer k,
	      float *a, Integer tda,   float *x, Integer incx );
extern void dtbsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n, Integer k,
	     double *a, Integer tda,  double *x, Integer incx );
extern void ctbsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n, Integer k,
	    Complex *a, Integer tda, Complex *x, Integer incx );
extern void ztbsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n, Integer k,
	    Zomplex *a, Integer tda, Zomplex *x, Integer incx );


/*               Solve a packed triangular system of equations               *
 *                                                                           *
 *                             y <- A**(-1)*x                                *
 *                             y <- A**(-T)*x        A is n-by-n triangular  *
 *                             y <- A**(-H)*x                in packed form  */

extern void stpsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	      float *ap,    float *x, Integer incx );
extern void dtpsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	     double *ap,   double *x, Integer incx );
extern void ctpsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	    Complex *ap,  Complex *x, Integer incx );
extern void ztpsv( MatrixTriangle uplo, MatrixTranspose trans,
	    MatrixUnitTriangular diag, Integer n,
	    Zomplex *ap,  Zomplex *x, Integer incx );

/*                               Rank-1 Updates                              *
 *                                                                           *
 *                            A <- alpha*x*yT + A          A is m-by-n real  */

extern void sger( Integer m, Integer n,     float alpha,
		   float *x, Integer incx,  float *y, Integer incy,
		   float *a, Integer tda );
extern void dger( Integer m, Integer n,    double alpha,
		  double *x, Integer incx, double *y, Integer incy,
		  double *a, Integer tda );


/*                            A <- alpha*x*yT + A         A is m-by-n complex */

extern void cgeru( Integer m, Integer n, Complex alpha, 
		  Complex *x, Integer incx,
	          Complex *y, Integer incy, Complex *a, Integer tda );
extern void zgeru( Integer m, Integer n, Zomplex alpha, 
		  Zomplex *x, Integer incx,
	          Zomplex *y, Integer incy, Zomplex *a, Integer tda );


/*                            A <- alpha*x*yH + A         A is m-by-n complex */

extern void cgerc( Integer m, Integer n, Complex alpha, 
		  Complex *x, Integer incx,
	          Complex *y, Integer incy, Complex *a, Integer tda );
extern void zgerc( Integer m, Integer n, Zomplex alpha, 
		  Zomplex *x, Integer incx,
	          Zomplex *y, Integer incy, Zomplex *a, Integer tda );

/*                            A <- alpha*x*xH + A  xx   A is n-by-n Hermitian */

extern void cher( MatrixTriangle uplo, Integer n, Complex alpha, 
	   Complex *x,   Integer incx, Complex *a, Integer tda );
extern void zher( MatrixTriangle uplo, Integer n, Zomplex alpha, 
	   Zomplex *x,   Integer incx, Zomplex *a, Integer tda );


/*                            A <- alpha*x*xH + A     A is n-by-n Hermitian   *
 *                                                           in packed form   */

extern void chpr( MatrixTriangle uplo, Integer n, Complex alpha, 
	   Complex *x,   Integer incx, Complex *ap );
extern void zhpr( MatrixTriangle uplo, Integer n, Zomplex alpha, 
	   Zomplex *x,   Integer incx, Zomplex *ap );


/*                           A <- alpha*x*x**T + A               A is n-by-n  *
								   symmetric  */

extern void ssyr( MatrixTriangle uplo, Integer n,  float alpha,
	     float *x, Integer incx,   float *a, Integer tda );
extern void dsyr( MatrixTriangle uplo, Integer n, double alpha,
	    double *x, Integer incx,  double *a, Integer tda );


/*                           A <- alpha*x*x**T + A                A is n-by-n *
							     symmetric packed */

extern void sspr( MatrixTriangle uplo, Integer n,  float alpha,
	     float *x, Integer incx,   float *ap );
extern void dspr( MatrixTriangle uplo, Integer n, double alpha,
	    double *x, Integer incx,  double *ap );


/*                               Rank-2 Updates                              *
 *                                                                           *
 *                   A <- alpha*x*y**H + y(alpha*x)**H + A       A is n-by-n *
 *							           Hermitian */

extern void cher2( MatrixTriangle uplo, Integer n, Complex alpha,
	   Complex *x,   Integer incx, Complex *y, Integer incy, 
	   Complex *a, Integer tda );
extern void zher2( MatrixTriangle uplo, Integer n, Zomplex alpha,
	   Zomplex *x,   Integer incx, Zomplex *y, Integer incy,
	   Zomplex *a, Integer tda );


/*                   A <- alpha*x*y**H + y(alpha*x)**H + A   A is n-by-n      *
							     Hermitian packed */

extern void chpr2( MatrixTriangle uplo, Integer n, Complex alpha,
	   Complex *x,   Integer incx, Complex *y, Integer incy, 
	   Complex *ap );
extern void zhpr2( MatrixTriangle uplo, Integer n, Zomplex alpha,
	   Zomplex *x,   Integer incx, Zomplex *y, Integer incy,
	   Zomplex *ap );


/*                   A <- alpha*x*y**T + alpha*y*x**T + A    A is n-by-n     *
 *								  symmetric  */

extern void ssyr2( MatrixTriangle uplo, Integer n,  float alpha,
	     float *x, Integer incx,   float *y, Integer incy,
	     float *a, Integer tda );
extern void dsyr2( MatrixTriangle uplo, Integer n, double alpha,
	    double *x, Integer incx,  double *y, Integer incy,
	    double *a, Integer tda );


/*                   A <- alpha*x*y**T + alpha*y*x**T + A    A is n-by-n      *
 *							     symmetric packed */

extern void sspr2( MatrixTriangle uplo, Integer n,  float alpha,
	     float *x, Integer incx,   float *y, Integer incy,
	     float *ap );
extern void dspr2( MatrixTriangle uplo, Integer n, double alpha,
	    double *x, Integer incx,  double *y, Integer incy,
	    double *ap );



/********                                                             ********
 ********                    LEVEL 3 BLAS                             ********
 ********                                                             ********/

/*                      Matrix - matrix multiplication                       *
 *                                                                           *
 *                     C <- alpha*op(A)*op(B) + beta*C                       *
 *   op(X) = X                                                  C is m-by-n  *
 *         = X**T                                                            *
 *         = X**H                                                            */

extern void sgemm( MatrixTranspose transa, MatrixTranspose transb,
	    Integer m, Integer n, Integer k,   float alpha,
	      float *a, Integer tda,   float *b, Integer tdb,
	      float beta,   float *c, Integer tdc );
extern void dgemm( MatrixTranspose transa, MatrixTranspose transb,
	    Integer m, Integer n, Integer k,  double alpha,
	     double *a, Integer tda,  double *b, Integer tdb,
	     double beta,  double *c, Integer tdc );
extern void cgemm( MatrixTranspose transa, MatrixTranspose transb,
	    Integer m, Integer n, Integer k, Complex alpha,
	    Complex *a, Integer tda, Complex *b, Integer tdb,
	    Complex beta, Complex *c, Integer tdc );
extern void zgemm( MatrixTranspose transa, MatrixTranspose transb,
	    Integer m, Integer n, Integer k, Zomplex alpha,
	    Zomplex *a, Integer tda, Zomplex *b, Integer tdb,
	    Zomplex beta, Zomplex *c, Integer tdc );


/*                 Symmetric Matrix - matrix multiplication                  *
 *                                                                           *
 *                          C <- alpha*A*B + beta*C     C is m-by-n          *
 *                          C <- alpha*B*A + beta*C     A is m-by-m symmetric*/

extern void ssymm( OperationSide side, MatrixTriangle uplo, 
		   Integer m,    Integer n,   float alpha,
		    float *a,    Integer tda, float *b, Integer tdb,
	            float beta,    float *c,  Integer tdc );
extern void dsymm( OperationSide side, MatrixTriangle uplo, 
		   Integer m,    Integer n,   double alpha,
		   double *a,    Integer tda, double *b, Integer tdb,
	           double beta,   double *c,  Integer tdc );
extern void csymm( OperationSide side, MatrixTriangle uplo, 
		   Integer m,    Integer n,   Complex alpha, 
		   Complex *a,   Integer tda, Complex *b, Integer tdb,
	           Complex beta, Complex *c,  Integer tdc );
extern void zsymm( OperationSide side, MatrixTriangle uplo, 
		   Integer m,    Integer n,   Zomplex alpha, 
		   Zomplex *a,   Integer tda, Zomplex *b, Integer tdb, 
		   Zomplex beta, Zomplex *c,  Integer tdc );


/*                 Hermitian Matrix - matrix multiplication                  *
 *                                                                           *
 *                          C <- alpha*A*B + beta*C    C is m-by-n           *
 *                          C <- alpha*B*A + beta*C    A is m-by-m Hermitian */

extern void chemm( OperationSide side, MatrixTriangle uplo, 
		   Integer m,    Integer n,   Complex alpha, 
		   Complex *a,   Integer tda, Complex *b, Integer tdb,
	           Complex beta, Complex *c,  Integer tdc );
extern void zhemm( OperationSide side, MatrixTriangle uplo, 
		   Integer m,    Integer n,   Zomplex alpha, 
		   Zomplex *a,   Integer tda, Zomplex *b, Integer tdb,
	           Zomplex beta, Zomplex *c, Integer tdc );


/*                               Rank-n Update                               *
 *                                                                           *
 *                        C <- alpha*A*A**T + beta*C                         *
 *                        C <- alpha*A**T*A + beta*C   C is n-by-n symmetric */

extern void ssyrk( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,      Integer k,
	            float alpha,     float *a, Integer tda,
	            float beta,      float *c, Integer tdc );
extern void dsyrk( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,      Integer k,
	            double alpha,   double *a, Integer tda,
	            double beta,    double *c, Integer tdc );
extern void csyrk( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,     Integer k,
	           Complex alpha, Complex *a, Integer tda,
	           Complex beta,  Complex *c, Integer tdc );
extern void zsyrk( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,     Integer k,
	           Zomplex alpha, Zomplex *a, Integer tda,
	           Zomplex beta,  Zomplex *c, Integer tdc );


/*                               Rank-n Update                               *
 *                                                                           *
 *                        C <- alpha*A*A**H + beta*C                         *
 *                        C <- alpha*A**H*A + beta*C   C is n-by-n Hermitian */

extern void cherk( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,     Integer k,
	           Complex alpha, Complex *a, Integer tda,
	           Complex beta,  Complex *c, Integer tdc );
extern void zherk( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,     Integer k,
	           Zomplex alpha, Zomplex *a, Integer tda,
	           Zomplex beta,  Zomplex *c, Integer tdc );


/*                 C <- alpha*A*B**H + alpha*B*A**H + beta*C     C is n-by-n *
 *                 C <- alpha*A**H*B + alpha*B**H*A + beta*C       Hermitian */

extern void cher2k( MatrixTriangle uplo, MatrixTranspose trans, 
		   Integer n,   Integer k,  Complex alpha, Complex *a, 
		   Integer tda, Complex *b, Integer tdb,   Complex beta,
		   Complex *c,  Integer tdc );
extern void zher2k( MatrixTriangle uplo, MatrixTranspose trans, 
		    Integer n,   Integer k,  Zomplex alpha, Zomplex *a, 
		    Integer tda, Zomplex *b, Integer tdb,   Zomplex beta,
		    Zomplex *c,  Integer tdc );


/*                 C <- alpha*A*B**T + alpha*B*A**T + beta*C     C is n-by-n *
 *                 C <- alpha*A**T*B + alpha*B**T*A + beta*C       symmetric */

extern void ssyr2k( MatrixTriangle uplo, MatrixTranspose trans, 
		    Integer n,   Integer k,    float alpha,  float *a, 
		    Integer tda,   float *b, Integer tdb,    float beta,
		      float *c,  Integer tdc );
extern void dsyr2k( MatrixTriangle uplo, MatrixTranspose trans, 
		    Integer n,   Integer k,   double alpha,  double *a, 
		    Integer tda,  double *b, Integer tdb,    double beta,
		     double *c,  Integer tdc );
extern void cher2k( MatrixTriangle uplo, MatrixTranspose trans, 
		    Integer n,   Integer k,  Complex alpha, Complex *a, 
		    Integer tda, Complex *b, Integer tdb,   Complex beta,
		    Complex *c,  Integer tdc );
extern void zher2k( MatrixTriangle uplo, MatrixTranspose trans, 
		    Integer n,   Integer k,  Zomplex alpha, Zomplex *a, 
		    Integer tda, Zomplex *b, Integer tdb,   Zomplex beta,
		    Zomplex *c,  Integer tdc );


/*               Triangular matrix - matrix multiplication                   *
 *                                                                           *
 *                             B <- alpha*op(A)B                             *
 *                             B <- alpha*B*op(A)                            *
 *      op(A) = A                                                            *
 *            = A**T                                  A is n-by-n triangular *
 *            = A**H                                  B is m-by-n            */

extern void strmm( MatrixTriangle side,    MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,    float alpha,   float *a, 
		   Integer tda,   float *b, Integer tdb );
extern void dtrmm( MatrixTriangle side,    MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,   double alpha,  double *a, 
		   Integer tda,  double *b, Integer tdb );
extern void ctrmm( MatrixTriangle side,    MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,  Complex alpha, Complex *a, 
		   Integer tda, Complex *b, Integer tdb );
extern void ztrmm( MatrixTriangle side,    MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,  Zomplex alpha, Zomplex *a, 
		   Integer tda, Zomplex *b, Integer tdb );


/*                  Solve a triangular system of equations                   *
 *                                                                           *
 *                        B <- alpha*op( A**(-1) )*B                         *
 *                        B <- alpha*B*op( A**(-1) )                         *
 *      op(A) = A                                                            *
 *            = A**T                                  A is n-by-n triangular *
 *            = A**H                                  B is m-by-n            */

extern void strsm( OperationSide side,     MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,    float alpha,   float *a, 
		   Integer tda,   float *b, Integer tdb );
extern void dtrsm( OperationSide side,     MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,   double alpha,  double *a, 
		   Integer tda,  double *b, Integer tdb );
extern void ctrsm( OperationSide side,     MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,  Complex alpha, Complex *a, 
		   Integer tda, Complex *b, Integer tdb );
extern void ztrsm( OperationSide side,     MatrixTriangle uplo,
	           MatrixTranspose transa, MatrixUnitTriangular diag,
	           Integer m,   Integer n,  Zomplex alpha, Zomplex *a, 
		   Integer tda, Zomplex *b, Integer tdb );

#endif /* !__CBLAS_H__ */

#ifdef __cplusplus
}
#endif

