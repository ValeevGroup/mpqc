/*!
   \file read.c
   \ingroup (PSIO)
*/
 
#include <unistd.h>
#include <string.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_READ(): Reads data from within a TOC entry from a PSI file.
**
**  \param unit   = The PSI unit number used to identify the file to all
**                  read and write functions.
**  \param key    = The TOC keyword identifying the desired entry.
**  \param buffer = The buffer to store the data as it is read.
**  \param size   = The number of bytes to read.
**  \param sadd   = The entry-relative starting page/offset of the desired data.
**  \param eadd   = A pointer to the entry-relative page/offset for the next
**                  byte after the end of the read request.
** 
** \ingroup (PSIO)
*/

int psio_read(unsigned int unit, char *key, char *buffer, ULI size,
	       psio_address sadd, psio_address *eadd)
{
  psio_ud *this_unit;
  psio_address address, end_address;
  psio_tocentry *this_entry;

  this_unit = &(psio_unit[unit]);
  
  /* Find the entry in the TOC */
  this_entry = psio_tocscan(unit, key);

  if(this_entry == NULL) {
      fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s\n", key);
      psio_error(unit,PSIO_ERROR_NOTOCENT);
  }
  else {

      /* Compute the starting page and offset for the block */
      address = psio_get_global_address(this_entry->sadd, sadd);

      /* Make sure the block starts and ends within the entry */
      if(address.page > this_entry->eadd.page)
	  psio_error(unit,PSIO_ERROR_BLKSTART);
      else if((address.page == this_entry->eadd.page) &&
	      (address.offset > this_entry->eadd.offset))
	  psio_error(unit,PSIO_ERROR_BLKSTART);

      end_address = psio_get_address(address, size);
      if((end_address.page > this_entry->eadd.page))
	  psio_error(unit,PSIO_ERROR_BLKEND);
      else if((end_address.page == this_entry->eadd.page) &&
	      (end_address.offset > this_entry->eadd.offset))
	  psio_error(unit,PSIO_ERROR_BLKEND);

      /* Update the eadd argument value for the caller */
      *eadd = psio_get_address(sadd, size);
    }

  /* Now read the actual data from the unit */
  psio_rw(unit, buffer, address, size, 0);

#ifdef PSIO_STATS
  psio_readlen[unit] += size;
#endif  

  return(0);
}

}
}
