/*!
** \file get_address.c
** \ingroup (PSIO)
*/
 
#include <util/psi3/libpsio/psio.h>

namespace psi3 {
namespace libpsio {

/*!
** PSIO_GET_ADDRESS(): Given a starting page/offset and a shift length
** (in bytes), return the page/offset of the next position in the file.
** \ingroup(PSIO)
*/

psio_address psio_get_address(psio_address start, ULI shift)
{
  psio_address address;
  ULI bytes_left;

  bytes_left = PSIO_PAGELEN - start.offset; /* Bytes remaining on fpage */

  if(shift >= bytes_left) { /* Shift to later page */
      address.page = start.page + (shift - bytes_left)/PSIO_PAGELEN + 1;
      address.offset = shift - bytes_left -
			(address.page - start.page - 1)*PSIO_PAGELEN;
    }
  else {  /* Block starts on current page */
      address.page = start.page;
      address.offset = start.offset + shift;
    }

  return address;
}

}
}
