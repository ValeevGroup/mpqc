
/* These routines manipulate keyword arrays.  A keyword
 * array differs from a data array in that its indices are indicated
 * by a keyword segment.  For example: array:0:1 = 6 is a keyword
 * array element.  The count is the upper bound plus 1. */

/* NOTE: If these routines are used to access keyword arrays, then
 * only the first place in the cwk in which a keyword array name is
 * found will be used. */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"
#include "ip_error.h"

#include "ip_error.gbl"

#include "ip_karray.gbl"
#include "ip_karray.lcl"

#include "ip_cwk.gbl"

GLOBAL_FUNCTION VOID
ip_cwk_karray_add_v(n,v)
int n;
int *v;
{
  int i;
  char indices[110],index[10];

  if (n>10) {
    ip_warn("ip_cwk_karray_add_v: too many indices: %d",n);
    return;
    }
  indices[0] = '\0';
  for (i=0; i<n; i++) {
    if (v[i] > 999999999 || v[i] < 0) {
      ip_warn("ip_cwk_karray_add_v: an index is too large or small: %d",v[i]);
      return;
      }
    sprintf(index,"%d",v[i]);
    strcpy(indices,index);
    if (i!=n-1) strcpy(indices,":");
    }
  ip_cwk_add(indices);
  return;
  }

GLOBAL_VA_FUNCTION VOID
#ifndef __STDC__
ip_cwk_karray_add(n)
int n;
#else
ip_cwk_karray_add(int n,...)
#endif
{
  va_list args;
  int i;
  int *v;

  if (n==0) {
    ip_cwk_karray_add_v(n,NULL);
    return;
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) {
      ip_warn("ip_cwk_karray_add: problem with malloc");
      return;
      }
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    free(v);
    ip_cwk_karray_add_v(n,v);
    return;
    }
  }

GLOBAL_FUNCTION ip_keyword_tree_t *
ip_karray_descend_v(kt,n,v)
ip_keyword_tree_t *kt;
int n;
int *v;
{
  ip_keyword_tree_t *r;
  int i;
  char index[10];

  if (!kt) return NULL;

  /* kt starts off at the array so we must first descend to the first
   * level of indices. */
  r = kt->down;
  if (!r) return NULL;

  for (i=0; i<n; i++) {
    if (!r) return r;
    if (v[i] > 999999999 || v[i] < 0) {
      ip_warn("ip_karray_descend_v: an index is too large or small: %d",v[i]);
      return NULL;
      }
    sprintf(index,"%d",v[i]);
    r = ip_descend_tree(r,index);
    if (r) r=r->down;
    }
  return r;
  }

GLOBAL_VA_FUNCTION ip_keyword_tree_t *
#ifndef __STDC__
ip_karray_descend(kt,n)
ip_keyword_tree_t *kt;
int n;
#else
ip_karray_descend(ip_keyword_tree_t *kt,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;

  if (n==0) {
    return ip_karray_descend_v(kt,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) {
      ip_warn("ip_karray_descend: problem with malloc");
      return NULL;
      }
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    free(v);
    return ip_karray_descend_v(kt,n,v);
    }
  }

GLOBAL_FUNCTION int
ip_count_v(keyword,karray_count,n,v)
char *keyword;
int *karray_count;
int n;
int *v;
{
  ip_keyword_tree_t *kt,*I;
  int index;
  int max;

  /* Descend the keyword tree to keyword. */
  kt = ip_cwk_descend_tree(keyword);
  if (kt == NULL) return IPE_KEY_NOT_FOUND;

  /* Descend the tree to the indices. */
  kt = ip_karray_descend_v(kt,n,v);
  if (kt == NULL) return IPE_NOT_AN_ARRAY;

  /* Go thru the keyword array and examine the indices. */
  I = kt;
  max = 0;
  do {
    if (sscanf(I->keyword,"%d",&index) != 1) return IPE_OUT_OF_BOUNDS;
    if (index<0) return IPE_OUT_OF_BOUNDS;
    if (index+1 > max) max = index + 1;
    } while ((I = I->across) != kt);

  *karray_count = max;
  return IPE_OK;
  }

/* This counts the number of elements in a keyword array. */
GLOBAL_VA_FUNCTION int
#ifndef __STDC__
ip_count(keyword,karray_count,n)
char *keyword;
int *karray_count;
int n;
#else
ip_count(char *keyword,int *karray_count,int n,...)
#endif
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_count_v(keyword,karray_count,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_count_v(keyword,karray_count,n,v);
    free(v);
    return r;
    }
  }
