/****************************************************************************
 * dirichlet.h
 * Author Doug Straub
 * Copyright 1991, Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Permission use, copy, and modify this software and its documentation
 * without fee for personal use or use within your organization is hereby
 * granted, provided that the above copyright notice is preserved in all
 * copies and that that copyright and this permission notice appear in
 * supporting documentation.  Permission to redistribute this software to
 * other organizations or individuals is not granted;  that must be
 * negotiated with the PSC.  Neither the PSC nor Carnegie Mellon
 * University make any representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *****************************************************************************/
/*
 * This module is a simple Gaussian elimination routine.
 */

#include "p3dgen.h"
#include "pgen_objects.h"
#include "dirichlet.h"
#include "ge_error.h"

static void row_div( float *array, int i, int j, int n, float f )
{
  int k;
  
  for ( k=j; k<=n; k++ ) 
    ( *( array+i*(n+1)+k ) ) /= f;
}

static void row_sub( float *array, int src1, int src2, int n, int start ) 
{
  int i;
  
  for ( i=start; i<=n; i++ ) 
    *( array+src2*( n+1 ) + i ) -= *( array+src1*( n+1 ) + i );
}


static void row_switch( float *array, int src1, int src2, int n ) 
{
  int i;
  float tmp;

  for ( i=0; i<=n; i++ ) {
    tmp = *( array+src2*( n+1 ) + i );
    *( array+src2*( n+1 ) + i ) = *( array+src1*( n+1 ) + i );
    *( array+src1*( n+1 ) + i ) = tmp;
  }
}

/*
 * Basic gaussian elimination.  if Ax=b, where A is an nxn matrix,
 * then array should be of the form: array = [ A b ].  The result is
 * returned through a pointer to an array, with a null pointer for
 * a non-invertible matrix.
 */ 
float *dch_gauss_elim( float *array, int n ) 
{
  int i, j;
  float *fptr;
  float *ans;

/* Forward elimination */

  for ( i=0; i<n; i++ ) {
    j = i;
    while ( ( ( *( array+ (n+1)*j + i ) ) == 0.0 ) && (j < n ) )
      j++;

    if ( j == n ) return( 0 );
    else if ( j != i )
      row_switch( array, i, j, n );

    row_div( array, i, i, n, *( array + (n+1)*i + i ) );
    for ( j=i+1; j<n; j++ ) { 
      fptr = array + (n+1)*j + i;
      if ( *fptr != 0.0 ) {
	row_div( array, j, i, n, *fptr );
	row_sub( array, i, j, n, i );  
      }
    }
  }
 
/* Back substitution */
 
  if ( !(ans = ( float * ) malloc( n*sizeof( float ) ) ) )
    ger_fatal("gauss_elim: unable to allocate %d floats!",n);

  j = 0;
  for ( i=( n-1 ); i>=0; i-- ) {
    *( ans + i ) = *( array + i*(n+1) + n);
    j= n - 1 - i;
    while ( j>0 ) {
      *( ans + i ) -= ( ( *( array + i*(n+1) + n-j) )*
			 ( *( ans + n-j ) ) );
      j--;
    }
  }
  return( ans );
}   

