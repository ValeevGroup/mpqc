
/* Storage integrals to avoid recomputation in iterative procedures.
 */

/* $Log$
 * Revision 1.3  1994/08/26 22:45:57  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:33:10  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:05:22  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:24:46  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.2  1992/01/30  01:31:36  cljanss
 * permutation information is now stored with the integrals
 *
 * Revision 1.1  1992/01/10  17:55:48  cljanss
 * Initial revision
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include "atoms.h"
#include "int_macros.h"
#include "int_flags.h"
#include "int_types.h"

#include "inter.h"

#define HASH_SIZE 4
#define MASK 0x3

struct intlist_struct {
  int sh1;
  int sh2;
  int sh3;
  int sh4;
  int size;
  int p12;
  int p34;
  int p13p24;
  double *buffer;
  struct intlist_struct *p;
  };
typedef struct intlist_struct intlist_t;

#include "storage.gbl"
#include "storage.lcl"

static intlist_t *hash_table[HASH_SIZE][HASH_SIZE][HASH_SIZE][HASH_SIZE];

static double *buffer,*available;

GLOBAL_FUNCTION VOID
int_storage(size)
int size;
{
  int i,j,k,l;
  if (!size) return;
  if (!int_maxsize) {
    fprintf(stderr,"int_maxsize is zero: call int_initialize first\n");
    fail();
    }
  buffer = (void *)0;
  while(size && !buffer) {
    available = buffer = (double *) malloc(sizeof(double)*size);
    if (!buffer) size /= 2;
    }
  if (!size) return;
  int_integral_storage = size;
  int_storage_threshold = int_maxsize / 4;
  int_used_integral_storage = 0;
  for (i=0; i<HASH_SIZE; i++) {
    for (j=0; j<HASH_SIZE; j++) {
      for (k=0; k<HASH_SIZE; k++) {
        for (l=0; l<HASH_SIZE; l++) {
          hash_table[i][j][k][l] = NULL;
          }
        }
      }
    }
  }

GLOBAL_FUNCTION VOID
int_reduce_storage_threshold()
{
  if (!int_integral_storage) return;
  int_storage_threshold = int_storage_threshold / 4;
  }

GLOBAL_FUNCTION VOID
int_done_storage()
{
  int i,j,k,l;
  if (!int_integral_storage) return;

  int_integral_storage = 0;
  free(buffer);
  buffer = NULL;
  available = NULL;
  int_used_integral_storage = 0;

  /* Free the allocated intlist storage. */
  for (i=0; i<HASH_SIZE; i++) {
    for (j=0; j<HASH_SIZE; j++) {
      for (k=0; k<HASH_SIZE; k++) {
        for (l=0; l<HASH_SIZE; l++) {
          free_intlist(hash_table[i][j][k][l]);
          }
        }
      }
    }
  }

GLOBAL_FUNCTION int
int_have_stored_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24)
int sh1;
int sh2;
int sh3;
int sh4;
int p12;
int p34;
int p13p24;
{
  int i;
  intlist_t *intlist = lookup(sh1,sh2,sh3,sh4,p12,p34,p13p24);

  if (!intlist) return 0;

#if 0
  printf("retrieving (%d %d|%d %d) size = %5d from 0x%x\n",
         sh1,sh2,sh3,sh4,intlist->size,intlist->buffer);
#endif
  for (i=0; i<intlist->size; i++) {
    int_buffer[i] = intlist->buffer[i];
    }
  
  return 1;
  }

GLOBAL_FUNCTION VOID
int_store_integral(sh1,sh2,sh3,sh4,p12,p34,p13p24,size)
int sh1;
int sh2;
int sh3;
int sh4;
int p12;
int p34;
int p13p24;
int size;
{

  if (lookup(sh1,sh2,sh3,sh4,p12,p34,p13p24)) {
    fprintf(stderr,"requested to store an integral already stored\n");
    fail();
    }
  if (int_integral_storage<size+int_used_integral_storage) return;
#if 0
  printf("storing (%d %d|%d %d) size = %5d at 0x%x\n",
         sh1,sh2,sh3,sh4,size,available);
#endif
  addint(sh1,sh2,sh3,sh4,p12,p34,p13p24,size);
  int_used_integral_storage += size;
  }

LOCAL_FUNCTION intlist_t *
lookup(sh1,sh2,sh3,sh4,p12,p34,p13p24)
int sh1;
int sh2;
int sh3;
int sh4;
int p12;
int p34;
int p13p24;
{
  intlist_t *intlist = hash_table[sh1&MASK][sh2&MASK][sh3&MASK][sh4&MASK];
  intlist_t *I;

  for (I=intlist; I!=NULL; I=I->p) {
    if (  (I->sh1 == sh1)
        &&(I->sh2 == sh2)
        &&(I->sh3 == sh3)
        &&(I->sh4 == sh4)
        &&(I->p12 == p12)
        &&(I->p34 == p34)
        &&(I->p13p24 == p13p24)) {
      return I;
      }
    }
  return NULL;
  }

LOCAL_FUNCTION VOID
addint(sh1,sh2,sh3,sh4,p12,p34,p13p24,size)
int sh1;
int sh2;
int sh3;
int sh4;
int p12;
int p34;
int p13p24;
int size;
{
  int i;
  intlist_t **intlist = &(hash_table[sh1&MASK][sh2&MASK][sh3&MASK][sh4&MASK]);
  intlist_t *tmp,*new;
  tmp = *intlist;
  new = (intlist_t *) malloc(sizeof(intlist_t));
  *intlist = new;
  new->sh1 = sh1;
  new->sh2 = sh2;
  new->sh3 = sh3;
  new->sh4 = sh4;
  new->p12 = p12;
  new->p34 = p34;
  new->p13p24 = p13p24;
  new->buffer = available;
  new->size = size; 
  new->p = tmp;
  for (i=0; i<size; i++) {
    *available++ = int_buffer[i];
    }
  }

LOCAL_FUNCTION VOID
free_intlist(intlist)
intlist_t *intlist;
{
  if (!intlist) return;
  free_intlist(intlist->p);
  free(intlist);
  }

LOCAL_FUNCTION VOID
fail()
{
  fprintf(stderr,"failing module:\n%s\n",__FILE__);
  exit(1);
  }

