/*!
   \file read_entry.c
   \ingroup (PSIO)
*/

#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_READ_ENTRY(): Reads an entire TOC entry from a PSI file.
**
**  \param unit   = The PSI unit number used to identify the file to all read
**                  and write functions.
**  \param key    = The TOC keyword identifying the desired entry.
**  \param buffer = The buffer to store the data as it is read.
**  \param size   = The number of bytes to read.
**
** Note that the value of size is not directly compared to the actual
** size of the entry, but care is taken to ensure that the end of the
** entry is not surpassed.
**
** \ingroup (PSIO)
*/

int psio_read_entry(unsigned int unit, char *key, char *buffer, ULI size)
{
  psio_address end;
  return psio_read(unit, key, buffer, size, PSIO_ZERO, &end);
}

}
}
