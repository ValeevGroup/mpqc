/*!
   \file write.c
   \ingroup (PSIO)
*/

#include <stdlib.h>
#include <string.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_WRITE(): Writes data to a TOC entry in a PSI file.
**
**  \param unit   = The PSI unit number used to identify the file to all read
**                  and write functions.
**  \param key    = The TOC keyword identifying the desired entry.
**  \param buffer = The buffer from which the data is written.
**  \param size   = The number of bytes to write.
**  \param sadd   = The entry-relative starting page/offset to write the data.
**  \param eadd   = A pointer to the entry-relative page/offset for the next
**                  byte after the end of the write request.
**
** \ingroup (PSIO)
*/

int psio_write(unsigned int unit, char *key, char *buffer, ULI size,
	       psio_address rel_start, psio_address *rel_end)
{
  psio_ud *this_unit;
  psio_tocentry *this_entry, *last_entry;
  psio_address address, end_address;

  this_unit = &(psio_unit[unit]);

  /* Find the entry in the TOC */
  this_entry = psio_tocscan(unit, key);

  if(this_entry == NULL) { /* New TOC entry */
    if(rel_start.page||rel_start.offset) psio_error(unit,PSIO_ERROR_BLKSTART);

    this_entry = (psio_tocentry *) malloc(sizeof(psio_tocentry));
    strcpy(this_entry->key,key);
    this_entry->next = NULL;
    this_entry->last = NULL;

    /* Compute the address of the entry */
    if(!(this_unit->toclen)) { /* First TOC entry */
      this_entry->sadd.page = 0;
      this_entry->sadd.offset = 3*sizeof(ULI);

      this_unit->toc = this_entry;
    }
    else {  /* Use ending address from last TOC entry */
      last_entry = psio_toclast(unit);
      this_entry->sadd = last_entry->eadd;

      last_entry->next = this_entry;
      this_entry->last = last_entry;
    }

    /* Data for the write call */
    address = this_entry->sadd;

    /* Set the end address for this_entry */
    this_entry->eadd = psio_get_address(this_entry->sadd, size);

    /* Update the unit's TOC stats */
    this_unit->toclen++;
    this_unit->tocaddress = this_entry->eadd;

    /* Update the rel_end argument value for the caller */
    *rel_end = psio_get_address(rel_start,size);
  }
  else { /* Old TOC entry */

    /* Compute the global starting page and offset for the block */
    address = psio_get_global_address(this_entry->sadd, rel_start);

    /* Make sure this block doesn't start past the end of the entry */
    if(address.page > this_entry->eadd.page)
      psio_error(unit,PSIO_ERROR_BLKSTART);
    else if((address.page == this_entry->eadd.page) &&
	    (address.offset > this_entry->eadd.offset))
      psio_error(unit,PSIO_ERROR_BLKSTART);

    /* Compute the new global ending address for the entry, if necessary */
    end_address = psio_get_address(address, size);
    if(end_address.page > this_entry->eadd.page) {
      if(this_entry->next != NULL) {
	fprintf(stderr, "PSIO_ERROR: Attempt to write into next entry: %d, %s\n", 
		unit, key);
	psio_error(unit, PSIO_ERROR_BLKEND);
      }
      this_entry->eadd = end_address;
      this_unit->tocaddress = end_address;
	  
    }
    else if((end_address.page == this_entry->eadd.page) &&
	    (end_address.offset > this_entry->eadd.offset))
      {
	if(this_entry->next != NULL) {
	  fprintf(stderr, "PSIO_ERROR: Attempt to write into next entry: %d, %s\n", 
		  unit, key);
	  psio_error(unit, PSIO_ERROR_BLKEND);
	}
	this_entry->eadd = end_address;
	this_unit->tocaddress = end_address;
      }

    /* Update the eadd argument value for the caller */
    *rel_end = psio_get_address(rel_start, size);
  }

  /* Now write the actual data to the unit */
  psio_rw(unit, buffer, address, size, 1);

#ifdef PSIO_STATS
  psio_writlen[unit] += size;
#endif

  return(0);
}

}
}
