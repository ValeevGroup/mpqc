#ifdef HAVE_CHEMISTRY_QC_CC
#  include <chemistry/qc/cc/LIBS.h>
#endif
#ifdef HAVE_CHEMISTRY_QC_PSI
#  include <chemistry/qc/psi/LIBS.h>
#endif
#ifdef HAVE_CHEMISTRY_QC_CINTS
#  include <chemistry/qc/cints/LIBS.h>
 #ifdef HAVE_CHEMISTRY_QC_MBPTR12
 #    include <chemistry/qc/mbptr12/LIBS.h>
 #endif
#endif
#ifdef HAVE_CHEMISTRY_CCA
#  include <chemistry/qc/intcca/LIBS.h>
#  include <util/misc/LIBS.h>
#endif
#include <chemistry/qc/mbpt/LIBS.h>
#include <chemistry/qc/dft/LIBS.h>
#include <chemistry/qc/scf/LIBS.h>
#include <util/options/LIBS.h>
