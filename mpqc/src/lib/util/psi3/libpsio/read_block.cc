/*!
   \file read_block.c
   \ingroup (PSIO)
*/

#include <unistd.h>
#include <string.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_READ_BLOCK(): Read a block of data.
**
** \param unit      = file number to read from
** \param key       = key to search for
** \param buffer    = where to put data
** \param blksiz    = ??
** \param start_blk = ??
** \param end_blk   = ??
**
** \ingroup (PSIO)
*/
int psio_read_block(unsigned int unit, char *key, char *buffer, ULI blksiz,
	            ULI start_blk, ULI end_blk)
{
  ULI size, shift;
  psio_ud *this_unit;
  psio_address sadd, eadd;
  psio_tocentry *this_entry;

  this_unit = &(psio_unit[unit]);
  
  /* Find the entry in the TOC */
  this_entry = psio_tocscan(unit, key);

  if(this_entry == NULL) {
    fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s\n", key);
    psio_error(unit,PSIO_ERROR_NOTOCENT);
  }
  else {
    size = (end_blk - start_blk + 1) * blksiz; /* The total buffer size */
    shift = start_blk * blksiz; /* Number of bytes to shift from start */

    /* Compute the starting page and offset for the block */
    sadd = psio_get_address(this_entry->sadd, shift);

    /* Make sure the block starts and ends within the entry */
    if((sadd.page > this_entry->eadd.page))
      psio_error(unit,PSIO_ERROR_BLKSTART);
    else if((sadd.page == this_entry->eadd.page) &&
	    (sadd.offset > this_entry->eadd.offset))
      psio_error(unit,PSIO_ERROR_BLKSTART);

    eadd = psio_get_address(sadd, size);
    if((eadd.page > this_entry->eadd.page))
      psio_error(unit,PSIO_ERROR_BLKEND);
    else if((eadd.page == this_entry->eadd.page) &&
	    (eadd.offset > this_entry->eadd.offset))
      psio_error(unit,PSIO_ERROR_BLKEND);

    /* Now read the actual data from the unit */
    psio_rw(unit, buffer, sadd, size, 0);
  }

  return(0);
}

}
}
