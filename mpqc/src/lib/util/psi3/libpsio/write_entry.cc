/*!
   \file write_entry.c
   \ingroup (PSIO)
*/

#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_WRITE_ENTRY()
**
** \ingroup (PSIO)
*/
int psio_write_entry(unsigned int unit, char *key, char *buffer, ULI size)
{
  psio_address end;
  return psio_write(unit, key, buffer, size, PSIO_ZERO, &end);
}

}
}
