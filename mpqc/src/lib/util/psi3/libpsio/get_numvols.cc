/*!
   \file get_numvols.c
   \ingroup (PSIO)
*/
 
#include <stdio.h>
#include <string.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_GET_NUMVOLS(): Get the number of volumes that file number 'unit'
** is split across.
**
** \ingroup (PSIO)
*/
unsigned int psio_get_numvols(unsigned int unit)
{
  unsigned int num;
  int errcod;
  char ip_token[PSIO_MAXSTR];

  num = 0;

  sprintf(ip_token,":PSI:FILES:FILE%u:NVOLUME",unit);
//  errcod = ip_data(ip_token,"%u",&num,0);
//  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":PSI:FILES:DEFAULT:NVOLUME");
//  errcod = ip_data(ip_token,"%u",&num,0);
//  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:NVOLUME",unit);
//  errcod = ip_data(ip_token,"%u",&num,0);
//  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NVOLUME");
//  errcod = ip_data(ip_token,"%u",&num,0);
//  if(errcod == IPE_OK) return(num);

  /* default to one volume */
  return(1);
}


/*!
** PSIO_GET_NUMVOLS_DEFAULT(): Get the number of volumes that file 
** number 'unit' is split across.
**
** \ingroup (PSIO)
*/
unsigned int psio_get_numvols_default(void)
{
  unsigned int num;
  int errcod;
  char ip_token[PSIO_MAXSTR];

  num = 0;

  sprintf(ip_token,":PSI:FILES:DEFAULT:NVOLUME");
//  errcod = ip_data(ip_token,"%u",&num,0);
//  if(errcod == IPE_OK) return(num);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:NVOLUME");
//  errcod = ip_data(ip_token,"%u",&num,0);
//  if(errcod == IPE_OK) return(num);

  /* default to one volume */
  return(1);
}

}
}
