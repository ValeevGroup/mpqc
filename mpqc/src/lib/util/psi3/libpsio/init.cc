/*!
   \file init.cc
   \ingroup (PSIO)
*/

#include <stdio.h>
#include <stdlib.h>

#include <util/psi3/libpsio/psio.h>
#include <util/psi3/libpsio/psifiles.h>

using namespace psi3::libpsio;

namespace psi3 {
namespace libpsio {

psio_ud *psio_unit;
psio_address PSIO_ZERO = {0,0};

/* Library state variable */
int _psi3_libpsio_state_ = 0;

#ifdef PSIO_STATS
ULI *psio_readlen;
ULI *psio_writlen;
#endif

/*!
** PSIO_INIT(): Allocates global memory needed by the I/O routines.
**
** No arguments.
**
** \ingroup (PSIO)
*/

int psio_init(void)
{
  int i,j;
  char *userhome;
  char filename[PSIO_MAXSTR];
  FILE *psirc;

  psio_unit = (psio_ud *) malloc(sizeof(psio_ud)*PSIO_MAXUNIT);

#ifdef PSIO_STATS
  psio_readlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
  psio_writlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
#endif

  if(psio_unit == NULL) {
    fprintf(stderr, "Error in PSIO_INIT()!\n");
    exit(PSI_RETURN_FAILURE);
  }

  for(i=0; i < PSIO_MAXUNIT; i++) {
#ifdef PSIO_STATS
    psio_readlen[i] = psio_writlen[i] = 0;
#endif      
    psio_unit[i].numvols = 0;
    for(j=0; j < PSIO_MAXVOL; j++) {
      psio_unit[i].vol[j].path = NULL;
      psio_unit[i].vol[j].stream = -1;
    }
    psio_unit[i].tocaddress.page = 0;
    psio_unit[i].tocaddress.offset = 0;
    psio_unit[i].toclen = 0;
    psio_unit[i].toc = NULL;
  }

  /* Open user's general .psirc file, if extant */
  userhome = getenv("HOME");
  sprintf(filename, "%s%s", userhome, "/.psirc");
  psirc = fopen(filename, "r");
  if(psirc != NULL) {
    // Need to parse the file elsewhere
    // ip_append(psirc, stdout);
    fclose(psirc);
  }

  /* Set library's state variable to initialized value (1) */
  _psi3_libpsio_state_ = 1;

  return(0);
}

/*!
** PSIO_STATE(): Returns state of the library (1=initialized, 0=noninitialized).
**
** No arguments.
**
** \ingroup (PSIO)
*/

int psio_state()
{
  return _psi3_libpsio_state_;
}

}
}

