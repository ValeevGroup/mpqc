
/* $Log$
 * Revision 1.3  1994/08/26 22:45:38  etseidl
 * fix a bunch of warnings, get rid of rcs id's, get rid of bread/bwrite and
 * fread/fwrite modules
 *
 * Revision 1.2  1993/12/30  13:32:55  etseidl
 * mostly rcs id stuff
 *
 * Revision 1.3  1992/06/17  22:04:50  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/31  01:22:44  jannsen
 * Merged in Sandia non CVS codes.
 *
 * Revision 1.1  1991/12/14  00:19:36  cljanss
 * Initial revision
 * */
/*
 * Revision 1.5  91/11/20  18:11:10  cljanss
 * Fixed an uninitialized variable bug.
 * 
 * Revision 1.4  91/10/31  14:46:42  cljanss
 * added the int_geometryinfo routine
 * 
 * Revision 1.3  91/09/28  19:26:52  cljanss
 * Switch to new file naming scheme
 * 
 * Revision 1.2  91/09/10  19:33:37  cljanss
 * int_basisinfo routine added.
 * 
 * Revision 1.1  1991/06/16  16:40:07  janssen
 * Initial revision
 * */

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include "atoms.h"
#include "int_macros.h"

#include "atomsprnt.h"
#include "inter.h"

#include "int_print.gbl"
#include "int_print.lcl"

/* Prints out information about the basis set. */
GLOBAL_FUNCTION VOID
int_basisinfo(fp,centers)
FILE *fp;
centers_t *centers;
{
  char *shell_name = "spdfghijklmnoqrtuvwxyz";
  int gc,c1,s1,i,j,k,index,offshell;
  index = 0;
  offshell=0;
  FOR_SHELLS(centers,c1,s1)
    FOR_GCCART2(gc,i,j,k,&(centers->center[c1].basis.shell[s1]))
      fprintf(fp," %3d: %c(%d,%d,%d) c = %2d (%s) s = %3d gc = %3d (%s)\n",
              index,shell_name[i+j+k],i,j,k,c1,
              centers->center[c1].atom,offshell,gc,
              centers->center[c1].basis.name);
      index++;
      END_FOR_GCCART2
    offshell++;
    END_FOR_SHELLS
  }

/* Print out information about the geometry and basis set. */
GLOBAL_FUNCTION VOID
int_geometryinfo(fp,centers)
FILE *fp;
centers_t *centers;
{
  int i,longest=0;
  char format[256];
  char *doubleformat;

  /* Find the longest basis name string. */
  for (i=0; i<centers->n; i++)
    if (longest<strlen(centers->center[i].basis.name))
      longest = strlen(centers->center[i].basis.name);

  /* Determine the print format. */
#if defined(NCUBE)
  doubleformat = "%14.8f";
#else
  doubleformat = "%14.8lf";
#endif
  sprintf(format," %%3d: %%2s %%%ds %s %s %s\n",
          longest,doubleformat,doubleformat,doubleformat);

  /* Print it out. */
  for (i=0; i<centers->n; i++) {
    fprintf(fp,format,
           i,centers->center[i].atom,
           centers->center[i].basis.name,
           centers->center[i].r[0],
           centers->center[i].r[1],
           centers->center[i].r[2]
           );
    }
  }


/* Prints out an integral buffer given
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * c1 = centers structure for center 1
 * s1 = shell number on center 1
 * ...
 * This prints out integrals using the offset arrays in the
 * centers structure.  Only nonzero integrals are printed.
 */
GLOBAL_FUNCTION VOID
int_offset_print(fp,buffer,c1,s1,c2,s2,c3,s3,c4,s4)
FILE *fp;
double *buffer;
centers_t *c1;
int s1;
centers_t *c2;
int s2;
centers_t *c3;
int s3;
centers_t *c4;
int s4;
{
  int nfunc1,nfunc2,nfunc3,nfunc4;

  nfunc1 = INT_SH_NFUNC(c1,s1);
  nfunc2 = INT_SH_NFUNC(c2,s2);
  nfunc3 = INT_SH_NFUNC(c3,s3);
  nfunc4 = INT_SH_NFUNC(c4,s4);

  int_offset_print_n(fp,buffer,nfunc1,nfunc2,nfunc3,nfunc4
    ,c1->func_offset + c1->func_num[s1]
    ,c2->func_offset + c2->func_num[s2]
    ,c3->func_offset + c3->func_num[s3]
    ,c4->func_offset + c4->func_num[s4]
    ,(c2==c1)&&(s2==s1)
    ,(c3==c1)&&(s3==s1) && (c4==c2)&&(s4==s2)
    ,(c4==c3)&&(s4==s3)
    );
  }

/* Prints out an integrals buffer given the number of functions
 * on each center and shell equivalency information.
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * n1 = number of functions in shell 1
 * ...
 * o1 = the basis function offset for shell1
 * ...
 * e12 = shell 1 == shell 2
 * e13e24 = (shell 1 == shell 3) && (shell 2 == shell 4)
 * e34 = shell 3 == shell 4
 */
GLOBAL_FUNCTION VOID
int_offset_print_n(fp,buffer,n1,n2,n3,n4,o1,o2,o3,o4,e12,e13e24,e34)
FILE *fp;
double *buffer;
int n1;
int n2;
int n3;
int n4;
int o1;
int o2;
int o3;
int o4;
int e12;
int e13e24;
int e34;
{
  int i,j,k,l;
  int index;

  index = 0;
  for (i=0; i<=INT_MAX1(n1); i++) {
    for (j=0; j<=INT_MAX2(e12,i,n2); j++) {
      for (k=0; k<=INT_MAX3(e13e24,i,n3); k++) {
        for (l=0; l<=INT_MAX4(e13e24,e34,i,j,k,n4); l++) {
          if (INT_NONZERO(buffer[index]))
            fprintf(fp," (%2d %2d|%2d %2d) = %11.7f\n",
                    o1+i,o2+j,o3+k,o4+l,buffer[index]);
          index++;
          }
        }
      }
    }
  }

/* Prints out an integral buffer given
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * c1 = centers structure for center 1
 * s1 = shell number on center 1
 * ...
 */
GLOBAL_FUNCTION VOID
int_print(fp,buffer,c1,s1,c2,s2,c3,s3,c4,s4)
FILE *fp;
double *buffer;
centers_t *c1;
int s1;
centers_t *c2;
int s2;
centers_t *c3;
int s3;
centers_t *c4;
int s4;
{
  int nfunc1,nfunc2,nfunc3,nfunc4;

  nfunc1 = INT_SH_NFUNC(c1,s1);
  nfunc2 = INT_SH_NFUNC(c2,s2);
  nfunc3 = INT_SH_NFUNC(c3,s3);
  nfunc4 = INT_SH_NFUNC(c4,s4);

  int_print_n(fp,buffer,nfunc1,nfunc2,nfunc3,nfunc4
    ,(c2==c1)&&(s2==s1)
    ,(c3==c1)&&(s3==s1) && (c4==c2)&&(s4==s2)
    ,(c4==c3)&&(s4==s3)
    );
  }

/* Prints out an integrals buffer given the number of functions
 * on each center and shell equivalency information.
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * n1 = number of functions in shell 1
 * ...
 * e12 = shell 1 == shell 2
 * e13e24  = (shell 1 == shell 3) && (shell 2 == shell 4)
 * e34 = shell 3 == shell 4
 */
GLOBAL_FUNCTION VOID
int_print_n(fp,buffer,n1,n2,n3,n4,e12,e13e24,e34)
FILE *fp;
double *buffer;
int n1;
int n2;
int n3;
int n4;
int e12;
int e13e24;
int e34;
{
  int i,j,k,l;
  int index;

  index = 0;
  for (i=0; i<=INT_MAX1(n1); i++) {
    for (j=0; j<=INT_MAX2(e12,i,n2); j++) {
      for (k=0; k<=INT_MAX3(e13e24,i,n3); k++) {
        for (l=0; l<=INT_MAX4(e13e24,e34,i,j,k,n4); l++) {
          if (INT_NONZERO(buffer[index]))
            fprintf(fp," (%2d %2d|%2d %2d) = (%4d) = %11.7f\n",
                    i,j,k,l,index,buffer[index]);
          index++;
          }
        }
      }
    }
  }

GLOBAL_FUNCTION VOID
int_print_intermediates(fp)
FILE *fp;
{
  fprintf(fp,"The integral intermediates:\n");

  fprintf(fp,"  int_prim_zeta:\n");
  print_double_matrix(fp,&int_prim_zeta);

  fprintf(fp,"  int_prim_k:\n");
  print_double_matrix(fp,&int_prim_k);

  fprintf(fp,"  int_prim_oo2zeta:\n");
  print_double_matrix(fp,&int_prim_oo2zeta);

  fprintf(fp,"  int_prim_p:\n");
  print_double_array3(fp,&int_prim_p);

  fprintf(fp,"  int_shell_to_prim:\n");
  print_int_vector(fp,&int_shell_to_prim);

  }
