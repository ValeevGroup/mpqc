/*!
   \file open.cc
   \ingroup (PSIO)
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_OPEN(): Opens a multivolume PSI direct access file for
** reading/writing data.
**
**  \param unit   = The PSI unit number used to identify the file to all
**                  read and write functions.
**  \param status = Indicates if the file is old (PSIO_OPEN_OLD) or new
**                  (PSIO_OPEN_NEW). 
**
** \ingroup (PSIO)
*/

int psio_open(unsigned int unit, int status)
{
  unsigned int i, j;
  int errcod, stream, tocstat;
  char name[PSIO_MAXSTR],path[PSIO_MAXSTR],fullpath[PSIO_MAXSTR];
  psio_ud *this_unit;

  this_unit = &(psio_unit[unit]);

  /* First check to see if this unit is aleady open */
  if(this_unit->vol[0].stream != -1) psio_error(unit,PSIO_ERROR_REOPEN);

  /* Get number of volumes to stripe across */
  this_unit->numvols = psio_get_numvols(unit);
  if(this_unit->numvols > PSIO_MAXVOL) psio_error(unit,PSIO_ERROR_MAXVOL);

  if(!(this_unit->numvols)) this_unit->numvols = 1;

  /* Get the file name prefix */
  errcod = psio_get_filename(unit,name);

  /* Build the name for each volume and open the file */
  for(i=0; i < this_unit->numvols; i++) {
    errcod = psio_get_volpath(unit, i, path);

    if(errcod && this_unit->numvols > 1)
      psio_error(unit,PSIO_ERROR_NOVOLPATH);

    sprintf(fullpath, "%s%s.%u", path, name, unit);
    this_unit->vol[i].path = (char *) malloc(strlen(fullpath)+1);
    strcpy(this_unit->vol[i].path,fullpath);

    /* Check if any previously opened volumes have the same path */
    for(j=0; j < i; j++)
      if (!strcmp(this_unit->vol[i].path,this_unit->vol[j].path))
	psio_error(unit,PSIO_ERROR_IDENTVOLPATH);

    /* Now open the volume */
    if(status == PSIO_OPEN_OLD) {
      this_unit->vol[i].stream =
	open(this_unit->vol[i].path,O_CREAT|O_RDWR,0644);
      if(this_unit->vol[i].stream == -1)
	psio_error(unit,PSIO_ERROR_OPEN);
    }
    else if(status == PSIO_OPEN_NEW) {
      this_unit->vol[i].stream =
	open(this_unit->vol[i].path,O_CREAT|O_RDWR|O_TRUNC,0644);
      if(this_unit->vol[i].stream == -1)
	psio_error(unit,PSIO_ERROR_OPEN);
    }
    else psio_error(unit,PSIO_ERROR_OSTAT);
  }

  if (status == PSIO_OPEN_OLD) tocstat = psio_tocread(unit);
  else if (status == PSIO_OPEN_NEW) {
      
    /* Init the TOC stats and write them to disk */
    this_unit->tocaddress.page = 0;
    this_unit->tocaddress.offset = 3*sizeof(ULI);
    this_unit->toclen = 0;
    this_unit->toc = NULL;

    /* Seek vol[0] to its beginning */
    stream = this_unit->vol[0].stream;
    errcod = lseek(stream, 0L, SEEK_SET);
    if(errcod == -1) psio_error(unit,PSIO_ERROR_LSEEK);
  
    errcod = write(stream, (char *) &(this_unit->tocaddress.page),
		   sizeof(ULI));
    if(errcod != sizeof(ULI)) psio_error(unit,PSIO_ERROR_WRITE);
    errcod = write(stream, (char *) &(this_unit->tocaddress.offset),
		   sizeof(ULI));
    if(errcod != sizeof(ULI)) psio_error(unit,PSIO_ERROR_WRITE);
    errcod = write(stream, (char *) &(this_unit->toclen), sizeof(ULI));
    if(errcod != sizeof(ULI)) psio_error(unit,PSIO_ERROR_WRITE);

  }
  else psio_error(unit,PSIO_ERROR_OSTAT);

  return(0);
}

}
}
