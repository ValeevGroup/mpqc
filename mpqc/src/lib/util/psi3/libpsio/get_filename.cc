/*!
   \file get_filename.c
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_GET_FILENAME(): Get the filename for filenumber 'unit'
**
** Returns: 
**   0 if a user-specified filename was found
**   1 if the global default will be used
**
** \ingroup (PSIO)
*/
int psio_get_filename(unsigned int unit, char *name)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":PSI:FILES:FILE%u:NAME",unit);
//  errcod = ip_data(ip_token,"%s",name,0);
//  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:NAME",unit);
//  errcod = ip_data(ip_token,"%s",name,0);
//  if(errcod == IPE_OK) return(0);

//  strcpy(name,psi_file_prefix);
  sprintf(name, "psiotest");
  return(1);
}


/*!
** PSIO_GET_FILENAME_DEFAULT(): Get the default filename
*/
int psio_get_filename_default(char *name)
{
//  strcpy(name,psi_file_prefix);
  sprintf(name, "psiotest");
  return(0);
}

}
}
