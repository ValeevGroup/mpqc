/*!
   \file tocwrite.c
   \ingroup (PSIO)
*/

#include <stdlib.h>
#include <unistd.h>
#include <util/psi3/libpsio/psio.h>
#include <util/psi3/libpsio/psifiles.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_TOCWRITE(): Write the table of contents for file number 'unit'.
**
** \ingroup (PSIO)
*/
int psio_tocwrite(unsigned int unit)
{
  unsigned int i;
  int errcod, volume, entry_size, stream;
  psio_ud *this_unit;
  psio_tocentry *this_entry;

  this_unit = &(psio_unit[unit]);
  entry_size = sizeof(psio_tocentry) - 2*sizeof(psio_tocentry *);

  /* Check that this unit is actually open */
  if(this_unit->vol[0].stream == -1) return(0);

  /* Seek vol[0] to its beginning */
  stream = this_unit->vol[0].stream;
  errcod = lseek(stream, 0L, SEEK_SET);
  if(errcod == -1) {
     fprintf(stderr, "Error in PSIO_TOCWRITE()!\n");
     exit(PSI_RETURN_FAILURE);
    }

  /* Dump the TOC stats to disk */
  errcod = write(stream, (char *) &(this_unit->tocaddress.page), sizeof(ULI));
  if(errcod != sizeof(ULI)) {
     fprintf(stderr, "Error in PSIO_TOCWRITE()!\n");
     exit(PSI_RETURN_FAILURE);
    }
  errcod = write(stream, (char *) &(this_unit->tocaddress.offset), sizeof(ULI));
  if(errcod != sizeof(ULI)) {
     fprintf(stderr, "Error in PSIO_TOCWRITE()!\n");
     exit(PSI_RETURN_FAILURE);
    }
  errcod = write(stream, (char *) &(this_unit->toclen), sizeof(ULI));
  if(errcod != sizeof(ULI)) {
     fprintf(stderr, "Error in PSIO_TOCWRITE()!\n");
     exit(PSI_RETURN_FAILURE);
    }

  /* Seek the TOC volume to the correct position */
  volume = (this_unit->tocaddress.page) % (this_unit->numvols);
  errcod = psio_volseek(&(this_unit->vol[volume]), this_unit->tocaddress.page,
	                this_unit->tocaddress.offset,this_unit->numvols);
  if(errcod == -1) {
      fprintf(stderr, "Error in PSIO_TOCWRITE()!\n");
      exit(PSI_RETURN_FAILURE);
    }

  this_entry = this_unit->toc;
  for(i=0; i < this_unit->toclen; i++) {

      /* This write() assumes a fixed ordering on the members of this_toc */
      errcod = write(this_unit->vol[volume].stream, (char *) this_entry,
		     entry_size);
      if(errcod != entry_size) {
      fprintf(stderr, "Error in PSIO_TOCWRITE()!\n");
      exit(PSI_RETURN_FAILURE);
        }

      this_entry = this_entry->next;
    }

  return(0);
}

}
}
