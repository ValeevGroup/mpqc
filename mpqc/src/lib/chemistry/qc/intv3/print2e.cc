
#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsprnt.h>
#include <chemistry/qc/intv3/macros.h>

#include <chemistry/qc/intv3/int2e.h>

/* Prints out an integral buffer given
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * c1 = centers structure for center 1
 * s1 = shell number on center 1
 * ...
 * This prints out integrals using the offset arrays in the
 * centers structure.  Only nonzero integrals are printed.
 */
void
Int2eV3::int_offset_print(FILE *fp,
                          double *buffer,
                          centers_t *c1, int s1,
                          centers_t *c2, int s2,
                          centers_t *c3, int s3,
                          centers_t *c4, int s4)
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
void
Int2eV3::int_offset_print_n(FILE *fp, double *buffer,
                            int n1, int n2, int n3, int n4,
                            int o1, int o2, int o3, int o4,
                            int e12, int e13e24, int e34)
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
void
Int2eV3::int_print(FILE *fp, double *buffer,
                   centers_t *c1, int s1,
                   centers_t *c2, int s2,
                   centers_t *c3, int s3,
                   centers_t *c4, int s4)
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
void
Int2eV3::int_print_n(FILE *fp, double *buffer,
                     int n1, int n2, int n3, int n4,
                     int e12, int e13e24, int e34)
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

void
Int2eV3::int_print_intermediates(FILE *fp)
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
