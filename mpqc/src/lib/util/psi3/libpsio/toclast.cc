/*!
   \file toclast.c
   \ingroup (PSIO)
*/

#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_TOCLAST(): Returns the last TOC entry.
**
** \ingroup (PSIO)
*/

psio_tocentry *psio_toclast(unsigned int unit)
{
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  while(this_entry->next != NULL) this_entry = this_entry->next;

  return(this_entry);
}

}
}
