libintv2.LIBSUF
#include <chemistry/molecule/LIBS.h>
#include <math/array/LIBS.h>
#include <util/keyval/LIBS.h>
#include <util/sgen/LIBS.h>
#include <chemistry/qc/intv2/MG.h>
#if MG != 0 && !defined(NO_TWO_ELECTRON)
#  include <chemistry/qc/oint2/LIBS.h>
#endif
