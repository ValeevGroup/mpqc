/*!
   \file tocclean.c
   \ingroup (PSIO)
*/

#include <string.h>
#include <stdlib.h>
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_TOCCLEAN(): Delete all TOC entries after the given key.
** If a blank key is given, the entire TOC will be wiped.
**
** \ingroup (PSIO)
*/

int psio_tocclean(unsigned int unit, char *key)
{
  psio_tocentry *this_entry, *last_entry, *prev_entry;

  /* Check the key length first */
  if((strlen(key)+1) > PSIO_KEYLEN) psio_error(unit,PSIO_ERROR_KEYLEN);

  this_entry = psio_tocscan(unit, key);
  if(this_entry == NULL) {
      if(!strcmp(key,"")) this_entry = psio_unit[unit].toc;
      else {
          fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s\n", key);
          psio_error(unit,PSIO_ERROR_NOTOCENT);
      }
    }
  else this_entry = this_entry->next;

  /* Get the end of the TOC and work backwards */
  last_entry = psio_toclast(unit);

  while((last_entry != this_entry) && (last_entry != NULL)) { 
      /* Now free all the remaining members */
      prev_entry = last_entry->last;
      free(last_entry);
      last_entry = prev_entry;
    }
  return(0);
}

}
}
