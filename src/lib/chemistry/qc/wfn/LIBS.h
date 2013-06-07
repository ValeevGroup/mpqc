libSCwfn.LIBSUF
#include <util/misc/LIBS.h>
#include <math/scmat/LIBS.h>
#include <math/mmisc/LIBS.h>
#include <chemistry/molecule/LIBS.h>
#include <chemistry/qc/basis/LIBS.h>
#include <chemistry/qc/intv3/LIBS.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS
#  include <chemistry/qc/cints/LIBS.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2
#  include <chemistry/qc/libint2/LIBS.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_CCA
#  include <chemistry/cca/int/LIBS.h>
#  include <chemistry/cca/server/LIBS.h>
#endif
