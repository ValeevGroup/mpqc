#include <bzerofast.h>

/* Commented out this version since the compiler
 * cannot handle it */
/*int bzerofast(double *d, int dimension)
  {
    int i;
  
    for (i=dimension; i; i--) {
      *d++ = 0.0;
      }
  
    return(0);
  }                     */
int bzerofast(double *d, int dimension)
{
  int i;

  for (i=0; i<dimension; i++) {
    *d++ = 0.0;
    }

  return(0);
}
