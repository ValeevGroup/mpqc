
#include <comm/picl/picl.h>

/* Gray code */
int gray0 (i)
     int i ;
{
  return( (i>>1)^i ) ;
}

/* inverse Gray code */
int ginv0 (i)
     int i ;
{
  int k ;
  k = i ;
  while ( k > 0 ) {
    k >>= 1 ;
    i ^= k ;
  }
  return (i) ;
}

