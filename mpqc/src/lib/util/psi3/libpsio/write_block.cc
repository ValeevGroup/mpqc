/*!
   \file write_block.c
   \ingroup (PSIO)
*/
 
#include <stdlib.h>
#include <string.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_WRITE_BLOCK(): 
**
** \ingroup (PSIO)
*/

int psio_write_block(unsigned int unit, char *key, char *buffer, ULI blksiz,
  		     ULI start_blk, ULI  end_blk)
{
  ULI size, shift;
  psio_ud *this_unit;
  psio_tocentry *this_entry, *last_entry;
  psio_address sadd, eadd;

  this_unit = &(psio_unit[unit]);

  /* Find the entry in the TOC */
  this_entry = psio_tocscan(unit, key);

  if(this_entry == NULL) { /* New TOC entry */
      if(start_blk) psio_error(unit,PSIO_ERROR_BLKSTART);

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
      size = (end_blk - start_blk + 1)*blksiz;
      sadd = this_entry->sadd;

      /* Set end address for this_entry */
      this_entry->eadd = psio_get_address(this_entry->sadd, size);

      /* Update the unit's TOC stats */
      this_unit->toclen++;
      this_unit->tocaddress = this_entry->eadd;
    }
  else { /* Old TOC entry */

      size = (end_blk - start_blk + 1) * blksiz; /* The total buffer size */
      shift = start_blk * blksiz; /* Number of bytes to shift from start */

      /* Compute the starting page and offset for the block */
      sadd = psio_get_address(this_entry->sadd, shift);

      /* Make sure this block doesn't start past the end of the entry */
      if(sadd.page > this_entry->eadd.page)
	  psio_error(unit,PSIO_ERROR_BLKSTART);
      else if((sadd.page == this_entry->eadd.page) &&
	      (sadd.offset > this_entry->eadd.offset))
	  psio_error(unit,PSIO_ERROR_BLKSTART);

      /* Compute the new ending address for the entry, if necessary */
      eadd = psio_get_address(sadd, size);
      if(eadd.page > this_entry->eadd.page) {
	  this_entry->eadd = this_unit->tocaddress = eadd;
	}
      else if((eadd.page == this_entry->eadd.page) &&
	      (eadd.offset > this_entry->eadd.offset))
	  {
	    this_entry->eadd = eadd;
	    this_unit->tocaddress = eadd;
	  }
    }

  /* Now write the actual data to the unit */
  psio_rw(unit, buffer, sadd, size, 1);

  return(0);
}

}
}
