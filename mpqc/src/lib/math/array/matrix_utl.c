
/* This file contains several useful routines for doing matrix operations. */

#include <stdio.h>
#include <tmpl.h>
#include "matrix.h"
#include "matrix_utl.gbl"
#include "matrix_utl.lcl"

/* Read a column ordered matrix on disk into a double matrix array.  The
 * storage for the matrix array must be preallocated.
 * fp is the file from which the matrix is read
 * matrix is the preallocated target matrix
 */
GLOBAL_FUNCTION void
read_col_double_matrix(fp,matrix)
FILE *fp;
double_matrix_t *matrix;
{
  int i,j;

  for (i=0; i<matrix->n2; i++) {
    for (j=0; j<matrix->n1; j++) {
      fread(&matrix->d[i][j],sizeof(double),1,fp);
      }
    }
  }

/* Fill the elements of a double matrix.
 * M is the input matrix
 * s1 = start row
 * e1 = end row + 1
 * s2 = start col
 * e2 = end col + 1
 * Mt is the preallocated output matrix
 * st1,et1,st2,et2 = the start and end rows and columns
 */
GLOBAL_FUNCTION void
fill_double_matrix(M,s1,e1,s2,e2, num)
double_matrix_t *M;
int s1;
int e1;
int s2;
int e2;
double num;
{
  int i,j;

  for (i=s1; i<e1; i++) {
    for (j=s2; j<e2; j++) {
      M->d[i][j] = num;
      }
    }
  }

/* Find the transpose of a double matrix.
 * M is the input matrix
 * s1 = start row
 * e1 = end row + 1
 * s2 = start col
 * e2 = end col + 1
 * Mt is the preallocated output matrix
 * st1,et1,st2,et2 = the start and end rows and columns
 */
GLOBAL_FUNCTION void
transpose_double_matrix(M,s1,e1,s2,e2, Mt,st1,et1,st2,et2)
double_matrix_t *M;
int s1;
int e1;
int s2;
int e2;
double_matrix_t *Mt;
int st1;
int et1;
int st2;
int et2;
{
  int i,j;
  int loop1,loop2;

  loop2 = st2;
  for (i=s1; i<e1; i++) {
    loop1 = st1;
    for (j=s2; j<e2; j++) {
      Mt->d[loop1][loop2] = M->d[i][j];
      loop1++;
      }
    loop2++;
    }
  }

/* This returns the sum of two matrices Apq += Bpq */
GLOBAL_FUNCTION VOID
add_double_matrix(A,sa1,ea1,sa2,ea2, B,sb1,eb1,sb2,eb2)
double_matrix_t *A;
int sa1;
int ea1;
int sa2;
int ea2;
double_matrix_t *B;
int sb1;
int eb1;
int sb2;
int eb2;
{
  int i,j,ib,jb;

  for (i=sa1,ib=sb1; i<ea1; i++,ib++) {
    for (j=sa2,jb=sb2; j<ea2; j++,jb++) {
      A->d[i][j] += B->d[ib][jb];
      }
    }
  }

/* This returns the dot product of two matrices d = Apq Bpq */
GLOBAL_FUNCTION double
dot_double_matrix(A,sa1,ea1,sa2,ea2, B,sb1,eb1,sb2,eb2)
double_matrix_t *A;
int sa1;
int ea1;
int sa2;
int ea2;
double_matrix_t *B;
int sb1;
int eb1;
int sb2;
int eb2;
{
  int i,j,ib,jb;
  double result = 0.0;

  for (i=sa1,ib=sb1; i<ea1; i++,ib++) {
    for (j=sa2,jb=sb2; j<ea2; j++,jb++) {
      result += A->d[i][j]*B->d[ib][jb];
      }
    }
  return result;
  }

/* This scales a matrix by a constant factor. */
GLOBAL_FUNCTION VOID
scale_double_matrix(A,factor)
double_matrix_t *A;
double factor;
{
  int i,j;
  for (i=0; i<A->n1; i++) {
    for (j=0; j<A->n2; j++) {
      A->d[i][j] *= factor;
      }
    }
  }

/* This computes the trace of a matrix. */
GLOBAL_FUNCTION double
trace_double_matrix(A)
double_matrix_t *A;
{
  double trace;
  int i;

  if (A->n1 != A->n2) {
    fprintf(stdout,"trace_matrix: tried to do trace of a nonsquare matrix\n");
    }

  trace = 0.0;
  for (i=0; i<A->n1; i++) {
    trace += A->d[i][i];
    }
  return trace;
  }

/* This multiplies two matrices.
 */
GLOBAL_FUNCTION void
multiply_double_matrix(A,sa1,ea1,sa2,ea2, B,sb1,eb1,sb2,eb2, C,sc1,ec1,sc2,ec2)
double_matrix_t *A;
int sa1;
int ea1;
int sa2;
int ea2;
double_matrix_t *B;
int sb1;
int eb1;
int sb2;
int eb2;
double_matrix_t *C;
int sc1;
int ec1;
int sc2;
int ec2;
{
  double tmp;
  int i,j,k;
  int ia,jb,kb;
  /* Check for consistency of the data. */
  if (  (ea1 > A->n1)
      ||(ea2 > A->n2)
      ||(eb1 > B->n1)
      ||(eb2 > B->n2)
      ||(ec1 > C->n1)
      ||(ec2 > C->n2)
      ||((ea2 - sa2) != (eb1 - sb1))
      ||((ea1 - sa1) != (ec1 - sc1))
      ||((eb2 - sb2) != (ec2 - sc2))) {
     printf("multiply_double_matrix: argument consistency error\n");
     printf(" A[%d:%d][%d:%d] * B[%d:%d][%d:%d] = C[%d:%d][%d:%d]\n",
	    sa1,ea1,sa2,ea2,sb1,eb1,sb2,eb2,sc1,ec1,sc2,ec2);
     printf(" real dimensions: A[%d][%d], B[%d][%d], C[%d][%d]\n",
	    A->n1,A->n2,B->n1,B->n2,C->n1,C->n2);
     exit(1);
     }

  ia = 0;
  for (i=sc1; i<ec1; i++) {
    jb = 0;
    for (j=sc2; j<ec2; j++) {
      tmp = 0.0;
      kb = 0;
      for (k=sa2; k<ea2; k++) {
        tmp += A->d[ia][k] * B->d[kb][jb];
        kb++;
        }
      C->d[i][j] = tmp;
      jb++;
      }
    ia++;
    }
  }
