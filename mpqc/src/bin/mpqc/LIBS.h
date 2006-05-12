#include <scdirlist.h>
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_CC
#  include <chemistry/qc/cc/LIBS.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI
#  include <chemistry/qc/psi/LIBS.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS
#  include <chemistry/qc/cints/LIBS.h>
 #ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12
 #    include <chemistry/qc/mbptr12/LIBS.h>
 #endif
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2
#  include <chemistry/qc/libint2/LIBS.h>
#endif
#ifdef HAVE_SC_SRC_LIB_CHEMISTRY_CCA
#  include <chemistry/cca/int/LIBS.h>
#  include <chemistry/cca/server/LIBS.h>
#  include <util/misc/LIBS.h>
#endif
#include <chemistry/qc/mbpt/LIBS.h>
#include <chemistry/qc/dft/LIBS.h>
#include <chemistry/qc/scf/LIBS.h>
#include <util/options/LIBS.h>
