/*!
   \file tocread.c
   \ingroup (PSIO)
*/

#include <unistd.h>
#include <stdlib.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_TOCREAD(): Read the table of contents for file number 'unit'.
** 
** \ingroup (PSIO)
*/
int psio_tocread(unsigned int unit)
{
  unsigned int i;
  int errcod, stream, volume, entry_size;
  psio_ud *this_unit;
  psio_tocentry *last_entry, *this_entry; 

  this_unit = &(psio_unit[unit]);
  entry_size = sizeof(psio_tocentry) - 2*sizeof(psio_tocentry *);

  /* Check that this unit is actually open */
  if(this_unit->vol[0].stream == -1) return(0);

  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  errcod = lseek(stream, 0L, SEEK_SET);
  if(errcod == -1) psio_error(unit,PSIO_ERROR_LSEEK);

  /* Read the TOC from disk */
  errcod = read(stream, (char *) &(this_unit->tocaddress.page),
                sizeof(ULI));
  if(errcod != sizeof(ULI)) return(1);
  errcod = read(stream, (char *) &(this_unit->tocaddress.offset),
                sizeof(ULI));
  if(errcod != sizeof(ULI)) return(1);
  errcod = read(stream, (char *) &(this_unit->toclen), sizeof(ULI));
  if(errcod != sizeof(ULI)) return(1);

  /* Malloc room for the TOC */
  this_unit->toc = (psio_tocentry *) malloc(sizeof(psio_tocentry));
  this_entry = this_unit->toc;
  this_entry->last = NULL;
  for(i=1; i < this_unit->toclen; i++) {
      last_entry = this_entry;
      this_entry = (psio_tocentry *) malloc(sizeof(psio_tocentry));
      last_entry->next = this_entry;
      this_entry->last = last_entry;
    }
  this_entry->next = NULL;

  /* Seek the TOC volume to the correct position */
  volume = (this_unit->tocaddress.page) % (this_unit->numvols);
  errcod = psio_volseek(&(this_unit->vol[volume]), this_unit->tocaddress.page,
	                this_unit->tocaddress.offset, this_unit->numvols);
  if(errcod == -1) psio_error(unit,PSIO_ERROR_LSEEK);

  /* Read the TOC entry-by-entry */
  this_entry = this_unit->toc;
  for(i=0; i < this_unit->toclen; i++) {

      /* This read() assumes a fixed ordering on the members of this_entry */
      errcod = read(this_unit->vol[volume].stream, (char *) this_entry,
		    entry_size);
      if(errcod != entry_size) psio_error(unit,PSIO_ERROR_READ);

      this_entry = this_entry->next;
    }

  return(0);
}

}
}
