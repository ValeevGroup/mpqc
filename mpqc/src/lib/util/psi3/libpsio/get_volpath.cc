/*!
   \file get_volpath.c
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <util/psi3/libpsio/psifiles.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*
** PSIO_GET_VOLPATH(): Get the path to a given volume for file number
** 'unit'.
**
** \ingroup (PSIO)
*/
int psio_get_volpath(unsigned int unit, unsigned int volume, char *path)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":PSI:FILES:FILE%u:VOLUME%u",unit,volume+1);
//  errcod = ip_data(ip_token,"%s",path,0);
//  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":PSI:FILES:DEFAULT:VOLUME%u",volume+1);
//  errcod = ip_data(ip_token,"%s",path,0);
//  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:FILE%u:VOLUME%u",unit,volume+1);
//  errcod = ip_data(ip_token,"%s",path,0);
//  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:VOLUME%u",volume+1);
//  errcod = ip_data(ip_token,"%s",path,0);
//  if(errcod == IPE_OK) return(0);

  /* default to /tmp/ for everything but chkpt */
  if(unit == PSIF_CHKPT) sprintf(path, "./");
  else sprintf(path, "/tmp/");
  return(1);
}


/*
** PSIO_GET_VOLPATH_DEFAULT(): Get the default path for the nth volume
** of any file.
**
** \ingroup (PSIO)
*/
int psio_get_volpath_default(unsigned int volume, char *path)
{
  int errcod;
  char ip_token[PSIO_MAXSTR];

  sprintf(ip_token,":PSI:FILES:DEFAULT:VOLUME%u",volume+1);
//  errcod = ip_data(ip_token,"%s",path,0);
//  if(errcod == IPE_OK) return(0);

  sprintf(ip_token,":DEFAULT:FILES:DEFAULT:VOLUME%u",volume+1);
//  errcod = ip_data(ip_token,"%s",path,0);
//  if(errcod == IPE_OK) return(0);

  /* default to /tmp/ */
  sprintf(path, "/tmp/");

  return(1);
}

}
}
